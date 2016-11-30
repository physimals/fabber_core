/*  inference_vb.h - VB inference technique class declarations

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once
#include "inference.h"
#include "convergence.h"

class PriorType
{
public:
	PriorType();
	PriorType(int idx, string param_name, FabberRunData &data);
	std::string m_param_name;
	int m_idx;
	char m_type;
	double m_prec;
	std::string m_filename;
	NEWMAT::RowVector m_image;

	void SetPrior(MVNDist *dist, int voxel);
};

std::ostream& operator<<(std::ostream& out, const PriorType value);

class VariationalBayesInferenceTechnique: public InferenceTechnique
{
public:
	/**
	 * Create a new instance of VariationalBayesInferenceTechnique
	 */
	static InferenceTechnique* NewInstance();

	virtual void GetOptions(std::vector<OptionSpec> &opts) const;
	virtual std::string GetDescription() const;
	virtual string GetVersion() const;

	VariationalBayesInferenceTechnique() :
			m_conv(NULL), initialFwdPrior(NULL), initialFwdPosterior(NULL), initialNoisePrior(NULL), initialNoisePosterior(
			NULL)
	{
	}

	virtual void Initialize(FwdModel* model, FabberRunData& args);
	virtual void DoCalculations(FabberRunData& data);
	virtual ~VariationalBayesInferenceTechnique();
protected:
	void InitializeMVNFromParam(FabberRunData& args, MVNDist *dist, string param_key);
	void InitializeNoiseFromParam(FabberRunData& args, NoiseParams *dist, string param_key);
	void MakeInitialDistributions(FabberRunData& args);
	void GetPriorTypes(FabberRunData& args);
	void LoadImagePriors(FabberRunData &allData);
	void PassModelData(int voxel);

	ConvergenceDetector* m_conv;

	// Initial priors and posteriors - not per voxel. Set up in Initialize
	std::auto_ptr<MVNDist> initialFwdPrior;
	MVNDist* initialFwdPosterior;
	NoiseParams* initialNoisePrior;
	NoiseParams* initialNoisePosterior;

	const NEWMAT::Matrix *m_origdata;
	const NEWMAT::Matrix *m_coords;
	const NEWMAT::Matrix *m_suppdata;
	int m_nvoxels;

	std::vector<PriorType> m_prior_types;

	// These are used for resuming a previous calculation
	std::string m_continueFromFile; // if empty, use initial posterior dists above
	std::string paramFilename;
	bool continueFwdOnly; // Only have fwd-model information

	// Reduce this to a linear problem, using the given
	// voxelwise linearizations (probably loaded from an MVN)
	std::string lockedLinearFile;

	bool m_printF;
	bool m_needF;
};

