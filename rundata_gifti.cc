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

