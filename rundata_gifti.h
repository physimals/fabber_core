/**
 * fabber_rundata_gifti.h
 *
 * Extends FabberRunData with loading and saving of data to
 * GIFTI files. This is the only part of Fabber which depends
 * on GIFTI and can be disabled if compiling for an environment
 * which has its own image load/save facilities.
 *
 * Martin Craig
 *
 * Copyright (C) 2019 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#ifndef NO_GIFTI
#include "rundata.h"

// There is a bug in the following header which means we MUST import fslsurfaceio.h before
// fslsurface.h (even though we don't need it here)
#include <fslsurface/fslsurfaceio.h>
#include <fslsurface/fslsurface.h>
#include <newmat.h>

#include <string>

/**
 * Run data which uses giftiio to load GIFTI files
 */
class FabberRunDataGifti : public FabberRunData
{
public:
    FabberRunDataGifti(bool compat_options = true);

    void Initialize();
    const NEWMAT::Matrix &LoadVoxelData(const std::string &filename);
    void SaveVoxelData(
        const std::string &filename, NEWMAT::Matrix &data, VoxelDataType data_type = VDT_SCALAR);
    void GetNeighbours(std::vector<std::vector<int> > &neighbours, 
                       std::vector<std::vector<int> > &neighbours2);
private:
    void SetCoordsFromSurface();
    fslsurface_name::fslSurface<float, unsigned int> m_surface;
};

#endif /* NO_NEWIMAGE */
