/*  inference_vb.h - VB inference technique class declarations

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once
#include "inference.h"

// Forward declaration -- see convergence.h, which is only actually
// #included in inference_vb.cc.  For now it's good enough just to know
// that the class exists. FIXME why not just include it?
class ConvergenceDetector;

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
		conv(NULL), initialFwdPrior(NULL), initialFwdPosterior(NULL), initialNoisePrior(NULL), initialNoisePosterior(
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

	ConvergenceDetector* conv;
	std::auto_ptr<MVNDist> initialFwdPrior;
	MVNDist* initialFwdPosterior;
	NoiseParams* initialNoisePrior;
	NoiseParams* initialNoisePosterior;

	const Matrix *m_origdata;
	const Matrix *m_coords;
	const Matrix *m_suppdata;

	// specification of priors from command line
	vector<RowVector> ImagePrior;
	vector<unsigned int> PSPidx;
	vector<char> PriorsTypes;
	vector<string> imagepriorstr;

	// These are used for resuming a previous calculation
	string continueFromFile; // if empty, use initial posterior dists above
	string paramFilename;
	bool continueFwdOnly; // Only have fwd-model information

	// Reduce this to a linear problem, using the given
	// voxelwise linearizations (probably loaded from an MVN)
	string lockedLinearFile;

	bool haltOnBadVoxel;
	bool printF;
	bool needF;
};

