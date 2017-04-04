/*  inference_vb.h - VB inference technique class declarations

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once
#include "convergence.h"
#include "inference.h"

class PriorType
{
public:
    PriorType();
    PriorType(int idx, std::vector<std::string> param_names, FabberRunData &data);
    std::string m_param_name;
    int m_idx;
    char m_type;
    double m_prec;
    std::string m_filename;
    NEWMAT::RowVector m_image;

    void SetPrior(MVNDist *dist, int voxel);

private:
    static std::string GetTypesString(FabberRunData &rundata, int num_params);
};

std::ostream &operator<<(std::ostream &out, const PriorType &value);

class VariationalBayesInferenceTechnique : public InferenceTechnique
{
public:
    /**
	 * Create a new instance of VariationalBayesInferenceTechnique
	 */
    static InferenceTechnique *NewInstance();

    virtual void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;
    virtual string GetVersion() const;

    VariationalBayesInferenceTechnique()
        : m_conv(NULL)
        , m_nvoxels(0)
        , m_noise_params(0)
        , m_origdata(NULL)
        , m_coords(NULL)
        , m_suppdata(NULL)
        , continueFwdOnly(false)
        , m_printF(false)
        , m_needF(false)
        , initialFwdPrior(NULL)
        , initialFwdPosterior(NULL)
        , initialNoisePrior(NULL)
        , initialNoisePosterior(NULL)
        , Nmcstep(0)
    {
    }

    virtual void Initialize(FwdModel *model, FabberRunData &args);
    virtual void DoCalculations(FabberRunData &data);
    virtual void SaveResults(FabberRunData &data) const;

    virtual ~VariationalBayesInferenceTechnique();

protected:
    void InitializeMVNFromParam(FabberRunData &args, MVNDist *dist, string param_key);
    void InitializeNoiseFromParam(FabberRunData &args, NoiseParams *dist, string param_key);
    void MakeInitialDistributions(FabberRunData &args);
    void GetPriorTypes(FabberRunData &args);
    void LoadImagePriors(FabberRunData &allData);
    void PassModelData(int voxel);

    /** Number of voxels in data */
    int m_nvoxels;

    /**
	 * Noise model in use. This is created by the inference
	 * method deleted on destruction
	 */
    std::auto_ptr<NoiseModel> noise;

    /**
	 * Number of noise parameters.
	 */
    int m_noise_params;

    /**
     * Convergence detector in use
     */
    ConvergenceDetector *m_conv;

    /** True if convergence detector requires the free energy */
    bool m_needF;

    /** True if we need to print the free energy at each iteration */
    bool m_printF;

    /** Prior types used for each model parameter */
    std::vector<PriorType> m_prior_types;

    // These are used for resuming a previous calculation
    /** If not empty, load initial MVN from this file */
    std::string m_continueFromFile;

    std::string paramFilename;

    /** If set, initial MVN only has fwd model information, not noise */
    bool continueFwdOnly;

    /** Free energy for each voxel */
    std::vector<double> resultFs;

    // Default priors and posteriors - not per voxel. Set up in Initialize
    // The fwd prior can be overridden by a per-voxel image prior
    // The fwd and noise posteriors can be overridden by a per-voxel MVN dist (restart run)
    std::auto_ptr<MVNDist> initialFwdPrior;
    MVNDist *initialFwdPosterior;
    NoiseParams *initialNoisePrior;
    NoiseParams *initialNoisePosterior;

    const NEWMAT::Matrix *m_origdata;
    const NEWMAT::Matrix *m_coords;
    const NEWMAT::Matrix *m_suppdata;

    /** Used by Adrian's spatial priors research */
    std::vector<MVNDist *> resultMVNsWithoutPrior;

    /**
     * Reduce this to a linear problem, using the given
     * voxelwise linearizations (probably loaded from an MVN)
     */
    std::string lockedLinearFile;

    /** Number of motion correction steps to run */
    int Nmcstep;
};
