/*  inference_spatialvb.h - implementation of VB with spatial priors

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "convergence.h"
#include "inference.h"
#include "run_context.h"

#include <string>
#include <vector>

class Vb : public InferenceTechnique
{
public:
    static InferenceTechnique *NewInstance();

    Vb()
        : m_nvoxels(0)
        , m_noise_params(0)
        , m_needF(false)
        , m_printF(false)
        , m_saveF(false)
        , m_origdata(NULL)
        , m_coords(NULL)
        , m_suppdata(NULL)
        , m_num_mcsteps(0)
        , m_spatial_dims(-1)
        , m_locked_linear(false)
    {
    }

    virtual void GetOptions(vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;
    virtual string GetVersion() const;

    virtual void Initialize(FwdModel *fwd_model, FabberRunData &args);
    virtual void DoCalculations(FabberRunData &data);

    virtual void SaveResults(FabberRunData &rundata) const;

protected:
    /**
     * Initialize noise prior or posterior distribution from a file stored in the
     * rundata under the given parameter key
     */
    void InitializeNoiseFromParam(FabberRunData &args, NoiseParams *dist, string param_key);

    /**
     * Pass the model the data, coords and suppdata for a voxel.
     *
     * FIXME this is not very nice and should not be necessary. Need to
     * audit what models are using this info and find alternatives, e.g.
     * reading suppdata in Initialize instead
     */
    void PassModelData(int voxel);

    /**
     * Determine whether we need spatial VB mode
     *
     * It is required either because it has been asked for (--method=spatialvb) or
     * if any spatial priors have been specified (types mMpP)
     */
    bool IsSpatial(FabberRunData &rundata) const;

    /**
     * Determine whether we are providing a Laplacian weighting matrix.
     */
    bool UseLaplacian(FabberRunData &rundata) const;

    /**
     * Do calculations loop in voxelwise mode (i.e. all iterations for
     * one voxel, then all iterations for the next voxel, etc)
     */
    virtual void DoCalculationsVoxelwise(FabberRunData &data);

    /**
     * Do calculations loop in spatial mode (i.e. one iteration of all
     * voxels, then next iteration of all voxels, etc)
     */
    virtual void DoCalculationsSpatial(FabberRunData &data);

    /**
     * Calculate free energy if required, and display if required
     */
    double CalculateF(int v, std::string label, double Fprior);

    /**
     * Output detailed debugging information for a voxel
     */
    void DebugVoxel(int v, const string &where);

    /**
     * Setup per-voxel data for Spatial VB
     *
     * Spatial VB needs each voxel's prior/posterior and other
     * data stored as it affects neighbouring voxels. This sets
     * up the vectors which store these things which are just
     * created on the fly for normal VB and throw away after each
     * voxel is done.
     */
    void SetupPerVoxelDists(FabberRunData &allData);

    /**
     * Ignore this voxel in future updates.
     *
     * No calculation of priors or posteriors will occur for this voxel
     * and it will be removed from the lists of neighbours for other voxels.
     * The effect should be as if it were masked
     */
    void IgnoreVoxel(int v);

    /** Number of voxels in data */
    int m_nvoxels;

    /**
     * Noise model in use. This is created by the inference
     * method deleted on destruction
     */
    std::auto_ptr<NoiseModel> m_noise;
    /**
     * Number of noise parameters.
     */
    int m_noise_params;

    /** True if convergence detector requires the free energy */
    bool m_needF;

    /** True if we need to print the free energy at each iteration */
    bool m_printF;

    /** True if we need to to save the final free energy */
    bool m_saveF;

    /** True if we need to to save the free energy history */
    bool m_saveFsHistory;

    /** Free energy for each voxel */
    std::vector<double> resultFs;

    /** Free energy history for each voxel / iteration number */
    std::vector<std::vector<double> > resultFsHistory;

    /** Voxelwise input data */
    const NEWMAT::Matrix *m_origdata;

    /** Voxelwise co-ordinates */
    const NEWMAT::Matrix *m_coords;

    /** Voxelwise supplementary data */
    const NEWMAT::Matrix *m_suppdata;

    /** NxN Laplacian weighting matrix */
    const NEWMAT::Matrix *m_laplacian;

    /** Number of motion correction steps to run */
    int m_num_mcsteps;

    /** Stores current run state (parameters, MVNs, linearization centres etc */
    RunContext *m_ctx;

    /** Linearized wrapper around the forward model */
    std::vector<LinearizedFwdModel> m_lin_model;

    /** Convergence detector for each voxel */
    std::vector<ConvergenceDetector *> m_conv;

    /**
     * Number of spatial dimensions
     *
     * 0 = no spatial smoothing
     * 1 = Probably not sensible!
     * 2 = Smoothing in slices only
     * 3 = Smoothing by volume
     */
    int m_spatial_dims;

    /**
     * Fix the linearization centres of the linearized forward model.
     *
     * This reduces the inference to a purely linear problem. The fixed
     * centres are generally loaded from an MVN file
     */
    bool m_locked_linear;
};
