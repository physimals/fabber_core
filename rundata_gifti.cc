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

#include <newmat.h>

#include <ostream>
#include <string>
#include <vector>

using namespace std;
using NEWMAT::Matrix;

FabberRunDataGifti::FabberRunDataGifti(bool compat_options)
    : FabberRunData(compat_options)
    , m_mask(1, 1, 1)
    , m_have_mask(false)
{
}

void FabberRunDataGifti::SetExtentFromSurface()
{
}

const Matrix &FabberRunDataGifti::LoadVoxelData(const std::string &filename)
{
    if (m_voxel_data.find(filename) == m_voxel_data.end())
    {
    }

    return m_voxel_data[filename];
}

void FabberRunDataGifti::SaveVoxelData(
    const std::string &filename, NEWMAT::Matrix &data, VoxelDataType data_type)
{
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

void FabberRunDataGifti::SetExtentFromSurface()
{
    LOG << "FabberRunDataGifti::Setting coordinates from surface" << endl;
}
