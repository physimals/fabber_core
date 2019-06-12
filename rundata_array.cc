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

void FabberRunDataArray::SetExtent(int nx, int ny, int nz, const int *mask)
{
    assert(nx > 0);
    assert(ny > 0);
    assert(nz > 0);
    assert(mask);

    m_extent.resize(3);
    m_extent[0] = nx;
    m_extent[1] = ny;
    m_extent[2] = nz;

    int nv = nx * ny * nz;
    if (mask)
    {
        m_mask.insert(m_mask.end(), mask, mask + nv);
    }
    else
    {
        m_mask.resize(nv, 1);
    }

    int *maskPtr = &m_mask[0];
    int v = 0;
    Matrix coords(3, nv);
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                bool masked = (*maskPtr == 0);
                if (!masked)
                {
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

void FabberRunDataArray::GetVoxelDataArray(string key, float *data)
{
    assert(data);
    Matrix mdata = FabberRunData::GetVoxelData(key);
    int nt = mdata.Nrows();
    float *dataPtr = data;

    for (int t = 0; t < nt; t++)
    {
        int v = 0;
        const int *maskPtr = &m_mask[0];
        for (int z = 0; z < m_extent[2]; z++)
        {
            for (int y = 0; y < m_extent[1]; y++)
            {
                for (int x = 0; x < m_extent[0]; x++)
                {
                    bool masked = (*maskPtr == 0);
                    if (!masked)
                        *dataPtr = mdata(t + 1, v + 1);
                    else
                        *dataPtr = 0;
                    ++dataPtr;
                    ++maskPtr;
                    if (!masked)
                        ++v;
                }
            }
        }
    }
}

void FabberRunDataArray::SetVoxelDataArray(string key, int data_size, const float *data)
{
    assert(data);
    const float *dataPtr = data;
    int num_voxels = m_extent[0] * m_extent[1] * m_extent[2];
    Matrix matrixData(data_size, num_voxels);

    int v = 0;
    for (int t = 0; t < data_size; t++)
    {
        v = 0;
        const int *maskPtr = &m_mask[0];
        for (int z = 0; z < m_extent[2]; z++)
        {
            for (int y = 0; y < m_extent[1]; y++)
            {
                for (int x = 0; x < m_extent[0]; x++)
                {
                    bool masked = (*maskPtr == 0);
                    if (!masked)
                    {
                        matrixData(t + 1, v + 1) = *dataPtr;
                    }
                    ++dataPtr;
                    ++maskPtr;
                    if (!masked)
                        ++v;
                }
            }
        }
    }
    matrixData = matrixData.Columns(1, v);
    FabberRunData::SetVoxelData(key, matrixData);
}

void FabberRunDataArray::GetNeighbours(std::vector<std::vector<int> > &neighbours, 
                                       std::vector<std::vector<int> > &neighbours2)
{
    LOG << "FabberRunDataArray::Getting nearest neigbours and second neighbours" << endl;
    int nx = m_extent[0];
    int ny = m_extent[1];
    int nz = m_extent[2];

    neighbours.clear();
    neighbours2.clear();

    // Tedious but fundamentally simple search for non-diagonal nearest neighbours
    int v = 0;
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
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
                                    (tz >= 0) && (tz < nz))
                                {
                                    // Note that voxel indices start at 1 because of NEWMAT
                                    neighbours[v].push_back(1+tx+ty*nx+tz*ny*nx);
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
}
