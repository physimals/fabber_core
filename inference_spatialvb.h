/*  inference_spatialvb.h - implementation of VB with spatial priors

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_vb.h"

class CovarianceCache : public Loggable
{
public:
    void CalcDistances(const NEWMAT::Matrix &voxelCoords, const string &distanceMeasure);
    const NEWMAT::SymmetricMatrix &GetDistances() const
    {
        return distances;
    }

    const NEWMAT::ReturnMatrix GetC(double delta) const; // quick to calculate
    const NEWMAT::SymmetricMatrix &GetCinv(double delta) const;

    //  const Matrix& GetCiCodist(double delta) const;
    const NEWMAT::SymmetricMatrix &GetCiCodistCi(double delta, double *CiCodistTrace = NULL) const;

    bool GetCachedInRange(double *guess, double lower, double upper, bool allowEndpoints = false) const;
    // If there's a cached value in (lower, upper), set *guess = value and
    // return true; otherwise return false and don't change *guess.

private:
    NEWMAT::SymmetricMatrix distances;
    typedef map<double, NEWMAT::SymmetricMatrix> Cinv_cache_type;
    mutable Cinv_cache_type Cinv_cache;
    mutable NEWMAT::SymmetricMatrix cinv;

    typedef map<double, pair<NEWMAT::SymmetricMatrix, double> > CiCodistCi_cache_type;
    //  mutable CiCodist_cache_type CiCodist_cache; // only really use the Trace
    mutable CiCodistCi_cache_type CiCodistCi_cache;
};

class SpatialVariationalBayes : public VariationalBayesInferenceTechnique
{
public:
    static InferenceTechnique *NewInstance();

    virtual void GetOptions(vector<OptionSpec> &opts) const;

    SpatialVariationalBayes()
        : VariationalBayesInferenceTechnique()
        , m_spatial_dims(-1)
        , m_spatial_speed(0)
        , m_shrinkage_type('-')
        , fixedDelta(0)
        , fixedRho(0)
        , m_update_first_iter(false)
        , m_use_evidence(false)
        , alwaysInitialDeltaGuess(false)
        , m_use_full_evidence(false)
        , m_use_sim_evidence(false)
        , m_alsoSaveWithoutPrior(false)
        , m_alsoSaveSpatialPriors(false)
        , m_lockedLinearEnabled(false)
        , firstParameterForFullEO(0)
        , useCovarianceMarginalsRatherThanPrecisions(false)
        , keepInterparameterCovariances(false)
        , newDeltaEvaluations(0)
        , bruteForceDeltaSearch(false)
    {
    }
    virtual void Initialize(FwdModel *fwd_model, FabberRunData &args);
    virtual void DoCalculations(FabberRunData &data);
    //    virtual ~SpatialVariationalBayes();

protected:
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
	 * Set up the StS matrix used for S and Z spatial priors
	 */
    void SetupStSMatrix();

    // Per-voxel prior and posterior distributions. For Spatial VB we need to
    // keep these around during iteration as the influence the calculations on
    // neighbouring voxels
    std::vector<NoiseParams *> noiseVox;      // these change. polymorphic type, so need to use pointers
    std::vector<NoiseParams *> noiseVoxPrior; // these may change in future
    std::vector<MVNDist> fwdPriorVox;
    std::vector<MVNDist> fwdPosteriorVox;
    std::vector<LinearizedFwdModel> linearVox;
    std::vector<MVNDist *> fwdPosteriorWithoutPrior;

    // StS matrix used for S and Z spatial priors
    NEWMAT::SymmetricMatrix StS;

    /**
	 * Number of spatial dimensions
	 *
	 * 0 = no spatial smoothing
	 * 1 = Probably not sensible!
	 * 2 = Smoothing in slices only
	 * 3 = Smoothing by volume
	 */
    int m_spatial_dims; // 0 = no spatial norm; 2 = slice only; 3 = volume

    //    bool useDataDrivenSmoothness;
    //    bool useShrinkageMethod;
    //    bool useDirichletBC;
    //    bool useMRF;
    //    bool useMRF2; // without the dirichlet bcs

    /**
	 * Maximum precision increase per iteration
	 */
    double m_spatial_speed; // Should be >1, or -1 = unlimited

    /**
	 * Type of spatial prior to use for each parameter. Should be one
	 * character per parameter, however if string ends with + then
	 * the last character is repeated for remaining parameters
	 */
    std::string m_prior_types_str;

    char m_shrinkage_type;

    /**
	 * Nearest-neighbours of each voxel. Vector size is number of voxels,
	 * each entry is vector of indices of nearest neighbours, starting at 1.
	 *
	 * FIXME Sparse matrix would be better?
	 */
    vector<vector<int> > m_neighbours;

    /**
	 * Next-nearest-neighbours of each voxel. Vector size is number of voxels,
	 * each entry is vector of indices of second nearest neighbours, starting at 1.
	 *
	 * FIXME Sparse matrix would be better?
	 */
    vector<vector<int> > m_neighbours2;

    /**
	 * Calculate first and second nearest neighbours of each voxel
	 */
    void CalcNeighbours(const NEWMAT::Matrix &voxelCoords);

    // For the new (Sahani-based) smoothing method:
    CovarianceCache covar;

    /**
	 * How to measure distances between voxels.
	 *
	 * dist1 = Euclidian distance
	 * dist2 = squared Euclidian distance
	 * mdist = Manhattan distance (|dx| + |dy|)
	 */
    std::string m_dist_measure;

    double fixedDelta;
    double fixedRho;

    /**
	 * Update spatial priors on first iteration?
	 */
    bool m_update_first_iter;

    /**
	 * Use evidence optimization
	 */
    bool m_use_evidence;
    double alwaysInitialDeltaGuess;

    /**
	 * Use full evidence optimization
	 */
    bool m_use_full_evidence;

    /**
	 * Use simultaneous evidence optimization
	 */
    bool m_use_sim_evidence;

    bool m_alsoSaveWithoutPrior;
    bool m_alsoSaveSpatialPriors;
    bool m_lockedLinearEnabled;

    int firstParameterForFullEO;
    bool useCovarianceMarginalsRatherThanPrecisions;
    bool keepInterparameterCovariances;

    int newDeltaEvaluations;

    bool bruteForceDeltaSearch;

    double OptimizeSmoothingScale(const NEWMAT::DiagonalMatrix &covRatio,
        //const SymmetricMatrix& covRatioSupplemented,
        const NEWMAT::ColumnVector &meanDiffRatio, double guess, double *optimizedRho = NULL, bool allowRhoToVary = true, bool allowDeltaToVary = true) const;

    double
    OptimizeEvidence(
        // const vector<MVNDist>& fwdPriorVox, // used for parameters other than k
        const vector<MVNDist *> &fwdPosteriorWithoutPrior, // used for parameter k
        // const vector<SymmetricMatrix>& Si,
        int k, const MVNDist *ifp, double guess, bool allowRhoToVary = false, double *rhoOut = NULL) const;
};
