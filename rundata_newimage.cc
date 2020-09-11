/**
 * fabber_rundata_newimage.cc
 *
 * Martin Craig
 *
 * Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "rundata_newimage.h"

#include "easylog.h"
#include "rundata.h"

#include <newimage/newimage.h>
#include <newimage/newimageio.h>
#include <newmat.h>

#include <ostream>
#include <string>
#include <vector>

using namespace std;
using namespace NEWIMAGE;
using NEWMAT::Matrix;

static void DumpVolumeInfo4D(const volume4D<float> &info, ostream &out)
{
    out << "FabberRunDataNewimage::Dimensions: x=" << info.xsize() << ", y=" << info.ysize()
        << ", z=" << info.zsize() << ", vols=" << info.tsize() << endl;
    out << "FabberRunDataNewimage::Voxel size: x=" << info.xdim() << "mm, y=" << info.ydim()
        << "mm, z=" << info.zdim() << "mm, TR=" << info.tdim() << " sec\n";
    out << "FabberRunDataNewimage::Intents: " << info.intent_code() << ", " << info.intent_param(1)
        << ", " << info.intent_param(2) << ", " << info.intent_param(3) << endl;
}

static void DumpVolumeInfo(const volume<float> &info, ostream &out)
{
    out << "FabberRunDataNewimage::Dimensions: x=" << info.xsize() << ", y=" << info.ysize()
        << ", z=" << info.zsize() << ", vols=1" << endl;
    out << "FabberRunDataNewimage::Voxel size: x=" << info.xdim() << "mm, y=" << info.ydim()
        << "mm, z=" << info.zdim() << "mm, TR=1"
        << " sec\n";
    out << "FabberRunDataNewimage::Intents: " << info.intent_code() << ", " << info.intent_param(1)
        << ", " << info.intent_param(2) << ", " << info.intent_param(3) << endl;
}

FabberRunDataNewimage::FabberRunDataNewimage(bool compat_options)
    : FabberRunData(compat_options)
    , m_mask(1, 1, 1)
    , m_have_mask(false)
{
}

void FabberRunDataNewimage::Initialize()
{
    string mask_fname = GetStringDefault("mask", "");
    m_have_mask = (mask_fname != "");

    if (m_have_mask)
    {
        LOG << "FabberRunDataNewimage::Initialize - Loading mask data from '" + mask_fname << "'" << endl;
        if (!fsl_imageexists(mask_fname))
        {
            throw DataNotFound(mask_fname, "File is invalid or does not exist");
        }
        read_volume(m_mask, mask_fname);
        m_mask.binarise(1e-16, m_mask.max() + 1, exclusive);
        DumpVolumeInfo(m_mask, LOG);
    }
        
    string data_fname = GetStringDefault("data", GetStringDefault("data1", ""));
    if (!fsl_imageexists(data_fname))
    {
        throw DataNotFound(data_fname, "File is invalid or does not exist");
    }
    read_volume(m_main_vol, data_fname);

    SetCoordsFromData();
}

const Matrix &FabberRunDataNewimage::LoadVoxelData(const std::string &filename)
{
    if (m_voxel_data.find(filename) == m_voxel_data.end())
    {
        // Load the data file using Newimage library
        if (!fsl_imageexists(filename))
        {
            throw DataNotFound(filename, "File is invalid or does not exist");
        }

        LOG << "FabberRunDataNewimage::Loading data from '" + filename << "'" << endl;
        volume4D<float> vol;
        try
        {
            read_volume4D(vol, filename);
            if (!m_have_mask)
            {
                // We need a mask volume so that when we save we can make sure
                // the image properties are set consistently with the source data
                m_mask = vol[0];
                m_mask = 1;
                m_have_mask = true;
            }
        }
        catch (...)
        {
            throw DataNotFound(filename, "Error loading file");
        }
        DumpVolumeInfo4D(vol, LOG);

        try
        {
            if (m_have_mask)
            {
                LOG << "FabberRunDataNewimage::Applying mask to data..." << endl;
                m_voxel_data[filename] = vol.matrix(m_mask);
            }
            else
            {
                m_voxel_data[filename] = vol.matrix();
            }
        }
        catch (exception &e)
        {
            LOG << "NEWMAT error while applying mask... Most likely a dimension mismatch. ***\n";
            throw;
        }
    }

    return m_voxel_data[filename];
}

void FabberRunDataNewimage::SaveVoxelData(
    const std::string &filename, NEWMAT::Matrix &data, VoxelDataType data_type)
{
    LOG << "FabberRunDataNewimage::Saving to nifti: " << filename << endl;
    int nifti_intent_code;
    switch (data_type)
    {
    case VDT_MVN:
        nifti_intent_code = NIFTI_INTENT_SYMMATRIX;
        break;
    default:
        nifti_intent_code = NIFTI_INTENT_NONE;
    }

    int data_size = data.Nrows();
    int nx = m_main_vol.xsize();
    int ny = m_main_vol.ysize();
    int nz = m_main_vol.zsize();
    volume4D<float> output(nx, ny, nz, data_size);
    if (m_have_mask)
    {
        output.setmatrix(data, m_mask);
    }
    else
    {
        output.setmatrix(data);
    }

    output.set_intent(nifti_intent_code, 0, 0, 0);
    output.setDisplayMaximumMinimum(output.max(), output.min());

    if (filename[0] == '/')
    {
        // Absolute path
        save_volume4D(output, filename);
    }
    else
    {
        // Relative path
        string filepath = GetOutputDir() + "/" + filename;
        save_volume4D(output, filepath);
    }
}

void FabberRunDataNewimage::SetCoordsFromData()
{
    LOG << "FabberRunDataNewimage::Setting coordinates from data" << endl;
    int nx = m_main_vol.xsize();
    int ny = m_main_vol.ysize();
    int nz = m_main_vol.zsize();

    volume4D<float> coordvol(nx, ny, nz, 3);
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                ColumnVector vcoord(3);
                vcoord << i << j << k;
                coordvol.setvoxelts(vcoord, i, j, k);
            }
        }
    }

    if (m_have_mask)
    {
        SetVoxelCoords(coordvol.matrix(m_mask));
    }
    else
    {
        SetVoxelCoords(coordvol.matrix());
    }
}

void FabberRunDataNewimage::GetNeighbours(std::vector<std::vector<int> > &neighbours, 
                                          std::vector<std::vector<int> > &neighbours2,
                                          std::vector<std::vector<double> > &weightings)
{
    LOG << "FabberRunDataNewimage::Getting nearest neigbours and second neighbours" << endl;
    neighbours.clear();
    neighbours2.clear();

    int nx = m_main_vol.xsize();
    int ny = m_main_vol.ysize();
    int nz = m_main_vol.zsize();

    // Mapping from xyz coords to masked voxel index
    vector<int> masked_idx(nx*ny*nz);
    int v = 0;
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                if (m_mask(x, y, z) != 0) 
                {
                    masked_idx[x+y*nx+z*ny*nx] = v;
                    v++;
                }
                else {
                    masked_idx[x+y*nx+z*ny*nx] = -1;
                }
            }
        }
    }

    // Tedious but fundamentally simple search for unmasked nearest neighbours
    v = 0;
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                if (m_mask(x, y, z) != 0) 
                {
                    neighbours.push_back(vector<int>());
                    for (int dx=-1; dx<2; dx++) 
                    {
                        for (int dy=-1; dy<2; dy++) 
                        {
                            for (int dz=-1; dz<2; dz++) 
                            {
                                // No diagonal neighbours
                                if (int(dx==0) + int(dy==0) + int(dz==0) == 2)
                                {
                                    int tx = x+dx;
                                    int ty = y+dy;
                                    int tz = z+dz;
                                    if ((tx >= 0) && (tx < nx) && 
                                        (ty >= 0) && (ty < ny) && 
                                        (tz >= 0) && (tz < nz) &&
                                        (m_mask(tx, ty, tz) != 0))
                                    {
                                        // Note that voxel indices start at 1 because of NEWMAT
                                        neighbours[v].push_back(1+masked_idx[tx+ty*nx+tz*ny*nx]);
                                    }
                                }
                            }
                        }
                    }
                    v++;
                }
            }
        }
    }
    int nv = v;
    neighbours2.resize(nv);

    // Search for Neighbours-of-neighbours, excluding self,
    // but including duplicates if there are two routes to get there (diagonally connected)
    // Note that voxel indices start at 1 because of NEWMAT
    for (unsigned int vid=1; vid <= nv; vid++)
    {
        // Go through the list of neighbours for each voxel.
        for (unsigned int n1 = 0; n1 < neighbours[vid-1].size(); n1++)
        {
            // n1id is the voxel index (starting at 1) of the neighbour
            unsigned int n1id = neighbours[vid-1][n1];
            int found_self_as_second_neighbour = 0;
            // Go through each of it's neighbours. Add each, apart from original voxel
            for (unsigned int n2=0; n2<neighbours[n1id-1].size(); n2++)
            {
                //cerr << n2 << ", " << n1id << ", " << neighbours.size() << ", " << nv << endl;
                unsigned int n2id = neighbours[n1id-1][n2];
                if (n2id != vid)
                {
                    neighbours2[vid-1].push_back(n2id);
                }
                else
                {
                    found_self_as_second_neighbour++;
                }
            }

            if (found_self_as_second_neighbour != 1)
            {
                throw FabberInternalError("Each voxel's neighbour must have "
                                          "the original voxel as a neighbour exactly once");
            }
        }
    }

    // For volumetric data, weightings are -1 for each neighbour and +<num nn> for the voxel self-weight
    weightings.clear();
    for (unsigned int vid=1; vid <= nv; vid++)
    {
        unsigned int nn = neighbours[vid-1].size();
        vector<double> voxel_weightings(nn+1, -1);
        voxel_weightings.back() = nn;
        weightings.push_back(voxel_weightings);
    }
}
