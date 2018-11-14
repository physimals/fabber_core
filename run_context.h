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

    /**
     * Start the processing
     *
     * This is a static function which calls the non-static `Run` method.
     * `obj` is a pointer to a subclass of ThreadContext. 
     *
     * This is required because thread creation (e.g. pthread_create) calls a 
     * function not a method. 
     */
    static void *Start(void *obj);

    FabberRunData *m_rundata;
    int m_id, m_num_workers;
    bool m_debug;
    bool m_halt_bad_voxel;
    NEWMAT::Matrix GetVoxelData(FabberRunData &rundata, std::string name);

    /**
     * Forward model to use.
     *
     * This is a per-thread copy of the model
     */
    FwdModel *m_model;

    /**
     * Noise model to use.
     *
     * This is a per-thread copy of the noise model
     */
    NoiseModel *m_noise;

    /** 
     * Voxelwise input data 
     *
     * This is a subset of the full input data from the start to the end
     * voxel
     */
    NEWMAT::Matrix m_origdata;

    /** 
     * Voxelwise co-ordinates 
     *
     * This is a subset of the coordinates, containing just those for the
     * subset of voxels this ThreadContext is processing
     */
    NEWMAT::Matrix m_coords;

    /**
     * Voxelwise supplementary data if provided
     *
     * This is a subset of the full suppdata for the voxels this ThreadContext
     * is processing. If no suppdata was provided it has 0 columns.
     */
    NEWMAT::Matrix m_suppdata;

    /**
     * Number of model parameters
     */
    int m_num_params;

    /**
     * Number of noise parameters
     */
    int m_noise_params;

    /** Current iteration, starting at 0 */
    int it;

    /** Current voxel, starting at 1 */
    int v;

    /** Start voxel number (in original data space) */
    int start_voxel;

    /** Total number of voxels to process in this ThreadContext */
    int nvoxels;

    bool m_locked_linear;

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

    /**
     * List of masked timepoints
     *
     * Masked timepoints are indexed starting at 1 and are ignored
     * in the analysis and parameter updates.
     */
    std::vector<int> m_masked_tpoints;

    void PassModelData(int v);
    void IgnoreVoxel(int v);
protected:

    /**
     * Run the calculation. 
     *
     * This must be implemented by subclasses
     */
    virtual void Run() {}

private:
    void InitializeNoiseFromParam(FabberRunData &rundata, NoiseParams *dist, std::string param_key);
    void InitMVNFromFile(FabberRunData &rundata);
};