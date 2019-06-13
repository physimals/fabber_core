/**
 * fabber_rundata_newimage.cc
 *
 * Martin Craig
 *
 * Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "rundata_gifti.h"

#include "easylog.h"
#include "rundata.h"

#include <fslsurface/fslsurfaceio.h>
#include <fslsurface/fslsurface.h>
#include <newmesh/newmesh.h>
#include <newmat.h>

#include <ostream>
#include <string>
#include <vector>

using namespace std;
using NEWMAT::Matrix;

FabberRunDataGifti::FabberRunDataGifti(bool compat_options)
    : FabberRunData(compat_options)
{
}

void FabberRunDataGifti::Initialize()
{
    string surface_fname = GetString("surface");
    LOG << "FabberRunDataGifti::Initialize - Reading surface file " << surface_fname << endl;
    fslsurface_name::read_surface(m_surface, surface_fname);
    LOG << "FabberRunDataGifti::Initialize - Vertices=" << m_surface.getNumberOfVertices() << ", faces=" << m_surface.getNumberOfFaces() << endl;
    SetCoordsFromSurface();
}

const Matrix &FabberRunDataGifti::LoadVoxelData(const std::string &filename)
{
    if (m_voxel_data.find(filename) == m_voxel_data.end())
    {
        LOG << "FabberRunDataGifti::LoadVoxelData filename '" + filename << "'" << endl;
        // Load the data file using fslsurface library
        // FIXME should check for presence of file before trying to load.
        //if (!fsl_imageexists(filename))
        //{
        //    throw DataNotFound(filename, "File is invalid or does not exist");
        //}

        fslsurface_name::fslSurface<float, unsigned int> vertex_data;
        try
        {
            fslsurface_name::read_surface(vertex_data, filename);
        }
        catch (...)
        {
            throw DataNotFound(filename, "Error loading surface data file");
        }

        LOG << "FabberRunDataGifti::LoadVoxelData - Scalars=" << vertex_data.getNumberOfScalarData() 
            << ", vectors=" << vertex_data.getNumberOfVectorData() 
            << endl;

        // We have no way of knowing which array the user wants to use so ensure there is only one
        vector<float> scalars = vertex_data.getScalars(0);
        if (vertex_data.getNumberOfScalarData() != 1) 
        {
            throw DataNotFound(filename, "Expected a single data array from Gifti file - found " + stringify(vertex_data.getNumberOfScalarData()));
        }

        // libfslsurface seems to write timeseries data as scalars with length NV*NT
        // We will use this convention for now but probably should raise this as a enhancement
        // to support timeseries data properly
        unsigned int nv = m_surface.getNumberOfVertices();
        if (scalars.size() % nv != 0) 
        {
            throw DataNotFound(filename, "Data not consistent with surface: " + stringify(scalars.size()) + 
                                         " values for " + stringify(nv) + " vertices");
        }
        unsigned int nt = scalars.size() / nv;
        NEWMAT::Matrix voxel_data(nt, nv);

        // We are making an assumption here about array ordering, but this is a result of 
        // libfslsurface munging a 2D array into 1D
        for (unsigned int v=0; v<nv; v++) 
        {
            for (unsigned int t=0; t<nt; t++) 
            {
                voxel_data(t+1, v+1) = scalars[v*nt+t];
            }
        }

        m_voxel_data[filename] = voxel_data;
    }

    return m_voxel_data[filename];
}

void FabberRunDataGifti::SaveVoxelData(
    const std::string &filename, NEWMAT::Matrix &data, VoxelDataType data_type)
{
    fslsurface_name::fslSurface<float, unsigned int> output_data;
    unsigned int nv = m_surface.getNumberOfVertices();
    if (data.Ncols() != nv) 
    {
        throw FabberInternalError("FabberRunDataGifti::Output data " + filename +
                                  " not consistent with reference surface: " + stringify(data.Ncols()) + 
                                  " values for " + stringify(nv) + " vertices");
    }
    unsigned int nt = data.Nrows();

    vector<float> vertex_data(nt*nv);
    for (unsigned int v=0; v<nv; v++) 
        {
            for (unsigned int t=0; t<nt; t++) 
            {
                vertex_data[v*nt+t] = data(t+1, v+1);
            }
        }

    output_data.addScalars(vertex_data, filename);

    if (filename[0] == '/')
    {
        // Absolute path
        fslsurface_name::writeGIFTI(output_data, filename + ".func.gii");
    }
    else
    {
        // Relative path
        string filepath = GetOutputDir() + "/" + filename;
        fslsurface_name::writeGIFTI(output_data, filepath + ".func.gii");
    }
}

void FabberRunDataGifti::SetCoordsFromSurface()
{
    LOG << "FabberRunDataGifti::Setting coordinates from surface" << endl;

    Matrix coords(3, m_surface.getNumberOfVertices());

    unsigned int idx = 1;
    for (std::vector< fslsurface_name::vertex<float> >::iterator iter=m_surface.vbegin();
         iter != m_surface.vend(); ++iter) 
    {
        coords(1, idx) = iter->x;
        coords(2, idx) = iter->y;
        coords(3, idx) = iter->z;
        ++idx;
    }

    SetVoxelCoords(coords);
}

void FabberRunDataGifti::GetNeighbours(std::vector<std::vector<int> > &neighbours, 
                                       std::vector<std::vector<int> > &neighbours2)
{
    LOG << "FabberRunDataGifti::Getting nearest neigbours and second neighbours" << endl;

    int nv = m_surface.getNumberOfVertices();
    int nt = m_surface.getNumberOfFaces();
    neighbours.clear();
    neighbours2.clear();
    neighbours.resize(nv);
    neighbours2.resize(nv);

    // Iterate over faces
    for (int t=0; t<nt; t++)
    {
        fslsurface_name::vec3<unsigned int> trig = m_surface.getFace(t, 3);

        // Iterate over vertices in a triangle - each vertex is connected to each other vertex
        // so add neighbours for all combinations unless they are already in the neighbours list
        // (no duplicates in the nearest neighbour lists)
        for (int v1=0; v1<3; v1++) 
        {
            if (std::find(neighbours[trig.x].begin(), neighbours[trig.x].end(), trig.y + 1) == neighbours[trig.x].end())
                neighbours[trig.x].push_back(trig.y + 1);
            if (std::find(neighbours[trig.x].begin(), neighbours[trig.x].end(), trig.z + 1) == neighbours[trig.x].end())
                neighbours[trig.x].push_back(trig.z + 1);
            if (std::find(neighbours[trig.y].begin(), neighbours[trig.y].end(), trig.x + 1) == neighbours[trig.y].end())
                neighbours[trig.y].push_back(trig.x + 1);
            if (std::find(neighbours[trig.y].begin(), neighbours[trig.y].end(), trig.z + 1) == neighbours[trig.y].end())
                neighbours[trig.y].push_back(trig.z + 1);
            if (std::find(neighbours[trig.z].begin(), neighbours[trig.z].end(), trig.x + 1) == neighbours[trig.z].end())
                neighbours[trig.z].push_back(trig.x + 1);
            if (std::find(neighbours[trig.z].begin(), neighbours[trig.z].end(), trig.y + 1) == neighbours[trig.z].end())
                neighbours[trig.z].push_back(trig.y + 1);
        }
    }

    // Search for Neighbours-of-neighbours, excluding self,
    // but including duplicates if there are two routes to get there (diagonally connected)
    // Note that voxel indices start at 1 because of NEWMAT
    for (unsigned int vid=1; vid <= nv; vid++)
    {   
        //cout << "Voxel " << (vid-1);
        // Go through the list of neighbours for each voxel.
        for (unsigned int n1 = 0; n1 < neighbours[vid-1].size(); n1++)
        {
            //cout << " " << (neighbours[vid-1][n1] - 1);
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
        //cout << endl;
    }
}

