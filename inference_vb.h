/*  inference_vb.h - VB inference technique class declarations

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once
#include "convergence.h"
#include "inference.h"

class PriorType {
public:
    PriorType();
    PriorType(int idx, string param_name, FabberRunData& data);
    std::string m_param_name;
    int m_idx;
    char m_type;
    double m_prec;
    std::string m_filename;
    NEWMAT::RowVector m_image;

    void SetPrior(MVNDist* dist, int voxel);
};

std::ostream& operator<<(std::ostream& out, const PriorType& value);

class VariationalBayesInferenceTechnique : public InferenceTechnique {
public:
    /**
	 * Create a new instance of VariationalBayesInferenceTechnique
	 */
    static InferenceTechnique* NewInstance();

    virtual void GetOptions(std::vector<OptionSpec>& opts) const;
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
        , m_outputOnly(false)
        , m_printF(false)
        , m_needF(false)
        , initialFwdPrior(NULL)
        , initialFwdPosterior(NULL)
        , initialNoisePrior(NULL)
        , initialNoisePosterior(
              NULL)
        , Nmcstep(0)
    {
    }

    virtual void Initialize(FwdModel* model, FabberRunData& args);
    virtual void DoCalculations(FabberRunData& data);
    virtual void SaveResults(FabberRunData& data) const;

    virtual ~VariationalBayesInferenceTechnique();

protected:
    void InitializeMVNFromParam(FabberRunData& args, MVNDist* dist, string param_key);
    void InitializeNoiseFromParam(FabberRunData& args, NoiseParams* dist, string param_key);
    void MakeInitialDistributions(FabberRunData& args);
    void GetPriorTypes(FabberRunData& args);
    void LoadImagePriors(FabberRunData& allData);
    void PassModelData(int voxel);

    /**
	 * Noise model in use. This is created by the inference
	 * method deleted on destruction
	 */
    std::auto_ptr<NoiseModel> noise;

    /**
	 * Number of noise parameters.
	 *
	 * This is used regularly so it's sensible to keep a
	 * copy around
	 */
    int m_noise_params;

    ConvergenceDetector* m_conv;

    /** Free energy for each voxel */
    std::vector<double> resultFs;

    // Initial priors and posteriors - not per voxel. Set up in Initialize
    std::auto_ptr<MVNDist> initialFwdPrior;
    MVNDist* initialFwdPosterior;
    NoiseParams* initialNoisePrior;
    NoiseParams* initialNoisePosterior;

    const NEWMAT::Matrix* m_origdata;
    const NEWMAT::Matrix* m_coords;
    const NEWMAT::Matrix* m_suppdata;
    int m_nvoxels;

    /**
	 * Used by Adrian's spatial priors research
	 */
    std::vector<MVNDist*> resultMVNsWithoutPrior;

    std::vector<PriorType> m_prior_types;

    // These are used for resuming a previous calculation
    std::string m_continueFromFile; // if empty, use initial posterior dists above
    std::string paramFilename;
    bool continueFwdOnly; // Only have fwd-model information
    bool m_outputOnly;

    // Reduce this to a linear problem, using the given
    // voxelwise linearizations (probably loaded from an MVN)
    std::string lockedLinearFile;

    bool m_printF;
    bool m_needF;

    /**
	 * Number of motion correction steps to run
	 */
    int Nmcstep;
};
