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
 * Class containing thread-specific state information for a run
 *
 * This includes the subset of data this particular thread is acting on,
 * a copy of the models being used (since Model at least is not thread safe
 * and the mutable output items used during inference (e.g. posteriors).
 *
 * The class also contains convenience accessors for voxel data which 
 * extract the subset pertinent to this thread.
 */
class ThreadContext : public Loggable
{
public:
    ThreadContext(FabberRunData &rundata, int worker_id=0, int n_workers=1, int start_vox=1, int num_vox=-1);
    ~ThreadContext();

    FabberRunData *m_rundata;
    int m_id, m_num_workers;
    bool m_debug;
    bool m_halt_bad_voxel;
    bool m_printF, m_needF, m_saveF;
    bool m_locked_linear;

    NEWMAT::Matrix GetVoxelData(FabberRunData &rundata, std::string name);

    FwdModel *m_model;
    NoiseModel *m_noise;

    /** Voxelwise input data */
    NEWMAT::Matrix m_origdata;

    /** Voxelwise co-ordinates */
    NEWMAT::Matrix m_coords;

    /** Voxelwise supplementary data */
    NEWMAT::Matrix m_suppdata;

    int m_num_params;
    int m_noise_params;

    /** Current iteration, starting at 0 */
    int it;

    /** Current voxel, starting at 1 */
    int v;

    /** Start voxel number (in original data space) */
    int start_voxel;

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

    void PassModelData(int v);
    void IgnoreVoxel(int v);
    /**
     * Calculate first and second nearest neighbours of each voxel
    */
    void CalcNeighbours(int spatial_dims);
private:
    void InitializeNoiseFromParam(FabberRunData &rundata, NoiseParams *dist, std::string param_key);
    void InitMVNFromFile(FabberRunData &rundata);
    
    /**
    * Check voxels are listed in order
    *
    * Order must be increasing in z value, or if same
    * increasing in y value, and if y and z are same
    * increasing in x value.
    *
    * This is basically column-major (Fortran) ordering - used as default by NEWIMAGE.
    */
    void CheckCoordMatrixCorrectlyOrdered();


};