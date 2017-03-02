/**
 * fabber_rundata_array.h
 *
 * Martin Craig
 *
 * Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "rundata.h"

#include <string>
#include <vector>

/**
 * Extends FabberRunData to allow setting and returning
 * data as plain C float arrays
 */
class FabberRunDataArray : public FabberRunData {
public:
    explicit FabberRunDataArray(bool compat_options = true)
        : FabberRunData(compat_options)
    {
    }

    void SetExtent(int nx, int ny, int nz, const int* mask);
    void GetVoxelData(std::string key, float* data);
    void SetVoxelData(std::string key, int data_size, const float* data);

private:
    std::vector<int> m_mask;
};