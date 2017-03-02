/**
 * fabber_rundata_array.cc
 *
 * Martin Craig
 *
 * Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "rundata_array.h"

#include "easylog.h"
#include "rundata.h"

#include "newmat.h"

#include <string>
#include <vector>

using namespace std;
using namespace NEWMAT;

void FabberRunDataArray::SetExtent(int nx, int ny, int nz, const int* mask)
{
    assert(nx > 0);
    assert(ny > 0);
    assert(nz > 0);
    assert(mask);

    FabberRunData::SetExtent(nx, ny, nz);

    int nv = nx * ny * nz;
    if (mask) {
        m_mask.insert(m_mask.end(), mask, mask + nv);
    } else {
        m_mask.resize(nv, 1);
    }

    int* maskPtr = &m_mask[0];
    int v = 0;
    Matrix coords(3, nv);
    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny; y++) {
            for (int z = 0; z < nz; z++) {
                bool masked = (*maskPtr == 0);
                if (!masked) {
                    coords(1, v + 1) = x;
                    coords(2, v + 1) = y;
                    coords(3, v + 1) = z;
                    ++v;
                }
                ++maskPtr;
            }
        }
    }

    coords = coords.Columns(1, v);
    SetVoxelCoords(coords);
}

void FabberRunDataArray::GetVoxelData(string key, float* data)
{
    assert(data);
    Matrix mdata = FabberRunData::GetVoxelData(key);
    int nt = mdata.Nrows();
    const int* maskPtr = &m_mask[0];
    float* dataPtr = data;

    int v = 0;
    for (int x = 0; x < m_extent[0]; x++) {
        for (int y = 0; y < m_extent[1]; y++) {
            for (int z = 0; z < m_extent[2]; z++) {
                bool masked = (*maskPtr == 0);
                for (int t = 0; t < nt; t++) {
                    if (!masked)
                        *dataPtr = mdata(t + 1, v + 1);
                    else
                        *dataPtr = 0;
                    ++dataPtr;
                }
                if (!masked)
                    ++v;
                ++maskPtr;
            }
        }
    }
}

void FabberRunDataArray::SetVoxelData(string key, int data_size, const float* data)
{
    assert(data);
    const int* maskPtr = &m_mask[0];
    const float* dataPtr = data;
    int num_voxels = m_extent[0] * m_extent[1] * m_extent[2];
    Matrix matrixData(data_size, num_voxels);

    int v = 0;
    for (int x = 0; x < m_extent[0]; x++) {
        for (int y = 0; y < m_extent[1]; y++) {
            for (int z = 0; z < m_extent[2]; z++) {
                bool masked = (*maskPtr == 0);
                for (int t = 0; t < data_size; t++) {
                    if (!masked) {
                        matrixData(t + 1, v + 1) = *dataPtr;
                    }
                    ++dataPtr;
                }
                if (!masked)
                    ++v;
                ++maskPtr;
            }
        }
    }
    matrixData = matrixData.Columns(1, v);
    FabberRunData::SetVoxelData(key, matrixData);
}
