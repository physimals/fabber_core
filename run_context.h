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
#include "convergence.h"

#include <vector>

/**
 * Structure containing mutable state information for the run
 *
 * This is currently just a convenient way to pass around the information needed.
 * However it could take on some of the roles of the VB inference technique, e.g.
 * storing a reference to the main voxel data and identifying nearest neighbours
 */
struct RunContext : public Loggable
{
    RunContext()
        : it(0)
        , v(1)
        , nvoxels(0)
    {
    }

    ~RunContext();

    void Initialize(FabberRunData &rundata);

    void PassModelData(int v);
    void InitializeNoiseFromParam(FabberRunData &rundata, NoiseParams *dist, std::string param_key);
    void InitMVNFromFile(FabberRunData &rundata, std::string paramFilename);
    
    /**
    * Check voxels are listed in order
    *
    * Order must be increasing in z value, or if same
    * increasing in y value, and if y and z are same
    * increasing in x value.
    *
    * This is basically column-major (Fortran) ordering - used as default by NEWIMAGE.
    */
    void CheckCoordMatrixCorrectlyOrdered(const NEWMAT::Matrix &voxelCoords);

    /**
     * Calculate first and second nearest neighbours of each voxel
    */
    void CalcNeighbours(const NEWMAT::Matrix &voxelCoords, int spatial_dims);

    FwdModel *m_model;
    NoiseModel *m_noise;

    /** Voxelwise input data */
    const NEWMAT::Matrix *m_origdata;

    /** Voxelwise co-ordinates */
    const NEWMAT::Matrix *m_coords;

    /** Voxelwise supplementary data */
    const NEWMAT::Matrix *m_suppdata;

    int m_num_params;
    int m_noise_params;

    /** Current iteration, starting at 0 */
    int it;

    /** Current voxel, starting at 1 */
    int v;

    /** Total number of voxels to process */
    int nvoxels;

    /** Voxels to ignore, indexed from 1 as per NEWMAT */
    std::vector<int> ignore_voxels;

    std::vector<MVNDist> fwd_prior;
    std::vector<MVNDist> fwd_post;
    std::vector<NoiseParams *> noise_prior;
    std::vector<NoiseParams *> noise_post;

    /** Free energy for each voxel */
    std::vector<double> resultFs;

    /**
     * Results of the inference method
     *
     * Vector of MVNDist, one for each voxel
     * Each MVNDist contains the means and covariance/precisions for
     * the parameters in the model
     */
    std::vector<MVNDist *> resultMVNs;

    /**
     * Linearized wrapper around the forward model - one for each voxel
     */
    std::vector<LinearizedFwdModel> m_lin_model;

    /**
     * Convergence detector for each voxel 
     */
    std::vector<ConvergenceDetector *> m_conv;

    std::vector<std::vector<int> > neighbours;
    std::vector<std::vector<int> > neighbours2;

    /**
     * List of masked timepoints
     *
     * Masked timepoints are indexed starting at 1 and are ignored
     * in the analysis and parameter updates.
     */
    std::vector<int> m_masked_tpoints;

};