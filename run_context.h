#pragma once
/**
 * run_context.h
 *
 * Copyright (C) 2007-2017 University of Oxford  
 */

/*  CCOPYRIGHT */

#include "noisemodel.h"
#include "fwdmodel_linear.h"
#include "dist_mvn.h"

#include <vector>

/**
 * Structure containing per-voxel state information
 * for the run
 */
struct RunContext
{
    RunContext(int nv) 
      : it(0), v(1), nvoxels(nv) {}
    
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