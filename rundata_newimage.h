/**
 * fabber_rundata_newimage.h
 *
 * Extends FabberRunData with loading and saving of data to
 * NIFTI files. This is the only part of Fabber which depends
 * on NEWIMAGE and can be disabled if compiling for an environment
 * which has its own image load/save facilities.
 *
 * Martin Craig
 *
 * Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#ifndef NO_NEWIMAGE
#include "rundata.h"

#include "newimage/newimage.h"
#include "newmat.h"

#include <string>
#include <vector>

/**
 * Run data which uses NEWIMAGE to load NIFTII files
 */
class FabberRunDataNewimage : public FabberRunData
{
public:
    FabberRunDataNewimage(bool compat_options = true);

    void Initialize();
    const NEWMAT::Matrix &LoadVoxelData(const std::string &filename);
    void SaveVoxelData(
        const std::string &filename, NEWMAT::Matrix &data, VoxelDataType data_type = VDT_SCALAR);
    void GetNeighbours(std::vector<std::vector<int> > &neighbours, 
                       std::vector<std::vector<int> > &neighbours2);

private:
    void SetCoordsFromData();
    NEWIMAGE::volume<float> m_mask;
    NEWIMAGE::volume<float> m_main_vol;
    bool m_have_mask;
};

#endif /* NO_NEWIMAGE */
