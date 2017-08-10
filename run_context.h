/**
 * run_context.h
 *
 * Copyright (C) 2007-2017 University of Oxford
 */

/*  CCOPYRIGHT */
#pragma once

#include "dist_mvn.h"
#include "fwdmodel_linear.h"
#include "noisemodel.h"

#include <vector>

/**
 * Structure containing per-voxel state information for the run
 *
 * This is currently just a convenient way to pass around the information needed.
 * However it could take on some of the roles of the VB inference technique, e.g.
 * storing a reference to the main voxel data and identifying nearest neighbours
 */
struct RunContext
{
    RunContext(int nv)
        : it(0)
        , v(1)
        , nvoxels(nv)
    {
    }

    /** Current iteration, starting at 0 */
    int it;

    /** Current voxel, starting at 1 */
    int v;

    /** Total number of voxels to process */
    int nvoxels;

    std::vector<MVNDist> fwd_prior;
    std::vector<MVNDist> fwd_post;
    std::vector<NoiseParams *> noise_prior;
    std::vector<NoiseParams *> noise_post;
    std::vector<std::vector<int> > neighbours;
    std::vector<std::vector<int> > neighbours2;
};