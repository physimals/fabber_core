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

#include "giftiio/giftiio.h"
#include "newmat.h"

#include <string>

/**
 * Run data which uses giftiio to load GIFTI files
 */
class FabberRunDataGifti : public FabberRunData
{
public:
    FabberRunDataGifti(bool compat_options = true);

    void SetExtentFromSurface();
    const NEWMAT::Matrix &LoadVoxelData(const std::string &filename);
    virtual void SaveVoxelData(
        const std::string &filename, NEWMAT::Matrix &data, VoxelDataType data_type = VDT_SCALAR);

private:
    void SetCoordsFromSurface(int nx, int ny, int nz);
};

#endif /* NO_NEWIMAGE */
