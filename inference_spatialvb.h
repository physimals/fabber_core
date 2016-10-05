/*  inference_spatialvb.h - implementation of VB with spatial priors

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_vb.h"

class CovarianceCache
{
public:

	void CalcDistances(const NEWMAT::Matrix& voxelCoords, const string& distanceMeasure);
	const NEWMAT::SymmetricMatrix& GetDistances() const
	{
		return distances;
	}

	const NEWMAT::ReturnMatrix GetC(double delta) const; // quick to calculate
	const NEWMAT::SymmetricMatrix& GetCinv(double delta) const;

	//  const Matrix& GetCiCodist(double delta) const;
	const NEWMAT::SymmetricMatrix& GetCiCodistCi(double delta, double* CiCodistTrace = NULL) const;

	bool GetCachedInRange(double *guess, double lower, double upper, bool allowEndpoints = false) const;
	// If there's a cached value in (lower, upper), set *guess = value and
	// return true; otherwise return false and don't change *guess.

private:
	NEWMAT::SymmetricMatrix distances;
	typedef map<double, NEWMAT::SymmetricMatrix> Cinv_cache_type;
	mutable Cinv_cache_type Cinv_cache;

	typedef map<double, pair<NEWMAT::SymmetricMatrix, double> > CiCodistCi_cache_type;
	//  mutable CiCodist_cache_type CiCodist_cache; // only really use the Trace
	mutable CiCodistCi_cache_type CiCodistCi_cache;
};

class SpatialVariationalBayes: public VariationalBayesInferenceTechnique
{
public:
	static InferenceTechnique* NewInstance();

	virtual void GetOptions(vector<OptionSpec> &opts) const;

	SpatialVariationalBayes() :
			VariationalBayesInferenceTechnique(), spatialDims(-1)
	{
	}
	virtual void Initialize(FwdModel* fwd_model, FabberRunData& args);
	virtual void DoCalculations(FabberRunData& data);
	//    virtual ~SpatialVariationalBayes();

protected:

	int spatialDims; // 0 = no spatial norm; 2 = slice only; 3 = volume
	bool continuingFromFile;

	//    bool useDataDrivenSmoothness;
	//    bool useShrinkageMethod;
	//    bool useDirichletBC;
	//    bool useMRF;
	//    bool useMRF2; // without the dirichlet bcs

	double maxPrecisionIncreasePerIteration; // Should be >1, or -1 = unlimited

	vector<vector<int> > neighbours; // Sparse matrix would be easier
	vector<vector<int> > neighbours2; // Sparse matrix would be easier
	void CalcNeighbours(const NEWMAT::Matrix& voxelCoords);

	//vector<string> imagepriorstr; now inherited from spatialvb

	// For the new (Sahani-based) smoothing method:
	CovarianceCache covar;
	string distanceMeasure;

	double fixedDelta;
	double fixedRho;
	bool updateSpatialPriorOnFirstIteration;
	bool useEvidenceOptimization;
	double alwaysInitialDeltaGuess;

	bool useFullEvidenceOptimization;
	bool useSimultaneousEvidenceOptimization;
	int firstParameterForFullEO;
	bool useCovarianceMarginalsRatherThanPrecisions;
	bool keepInterparameterCovariances;

	int newDeltaEvaluations;

	string spatialPriorsTypes; // one character per parameter
	//    bool spatialPriorOutputCorrection;

	bool bruteForceDeltaSearch;

	double OptimizeSmoothingScale(const NEWMAT::DiagonalMatrix& covRatio,
	//const SymmetricMatrix& covRatioSupplemented,
			const NEWMAT::ColumnVector& meanDiffRatio, double guess, double* optimizedRho = NULL, bool allowRhoToVary =
					true, bool allowDeltaToVary = true) const;

	double
	OptimizeEvidence(
	// const vector<MVNDist>& fwdPriorVox, // used for parameters other than k
			const vector<MVNDist*>& fwdPosteriorWithoutPrior, // used for parameter k
			// const vector<SymmetricMatrix>& Si,
			int k, const MVNDist* ifp, double guess, bool allowRhoToVary = false, double* rhoOut =
			NULL) const;
};

