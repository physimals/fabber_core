/*  inference_spatialvb.h - implementation of VB with spatial priors

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_vb.h"

#include "covariance_cache.h"

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
        , m_fixed_delta(0)
        , m_fixed_rho(0)
        , m_update_first_iter(false)
        , m_use_evidence(false)
        , m_always_inital_delta_guess(false)
        , m_use_full_evidence(false)
        , m_use_sim_evidence(false)
        , m_save_without_prior(false)
        , m_save_spatial_priors(false)
        , m_locked_linear(false)
        , m_full_eo_first_param(0)
        , m_use_covar_marginals_not_precisions(false)
        , m_keep_param_covars(false)
        , m_new_delta_evaluations(0)
        , m_brute_force_delta_search(false)
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
	 * Set up the StS matrix used for S and Z spatial priors
	 */
    void SetupStSMatrix();

    void UpdateAkmean(NEWMAT::DiagonalMatrix &akmean);
    void UpdateDeltaRho(NEWMAT::DiagonalMatrix &delta, NEWMAT::DiagonalMatrix &rho,
        const NEWMAT::DiagonalMatrix &akmean, bool first_iter);
    void CalculateCinv(std::vector<NEWMAT::SymmetricMatrix> &Sinvs, NEWMAT::DiagonalMatrix &delta,
        NEWMAT::DiagonalMatrix &rho, NEWMAT::DiagonalMatrix &akmean);
    void DoSimEvidence(std::vector<NEWMAT::SymmetricMatrix> &Sinvs);
    void DoFullEvidence(std::vector<NEWMAT::SymmetricMatrix> &Sinvs);
    void SetFwdPriorShrinkageTypeS(int voxel, const NEWMAT::DiagonalMatrix &akmean);
    void SetFwdPriorShrinkageType(int voxel, const NEWMAT::DiagonalMatrix &akmean);
    double SetFwdPrior(int voxel, const std::vector<NEWMAT::SymmetricMatrix> &Sinvs, bool isFirstIteration);

    /**
	 * Calculate first and second nearest neighbours of each voxel
	*/
    void CalcNeighbours(const NEWMAT::Matrix &voxelCoords);

    double OptimizeSmoothingScale(const NEWMAT::DiagonalMatrix &covRatio,
        const NEWMAT::ColumnVector &meanDiffRatio, double guess, double *optimizedRho = NULL, bool allowRhoToVary = true, bool allowDeltaToVary = true) const;

    double OptimizeEvidence(
        const std::vector<MVNDist *> &fwdPosteriorWithoutPrior, // used for parameter k
        int k, const MVNDist *ifp, double guess, bool allowRhoToVary = false, double *rhoOut = NULL) const;

    // Per-voxel prior and posterior distributions. For Spatial VB we need to
    // keep these around during iteration as the influence the calculations on
    // neighbouring voxels
    // NoiseParams is polymorphic type, so need to use pointers
    std::vector<NoiseParams *> m_noise_post;
    std::vector<NoiseParams *> m_noise_prior;
    std::vector<MVNDist> m_fwd_prior;
    std::vector<MVNDist> m_fwd_post;
    std::vector<LinearizedFwdModel> m_lin_model;
    std::vector<MVNDist *> m_fwd_post_no_prior;

    /**
	 * Number of spatial dimensions
	 *
	 * 0 = no spatial smoothing
	 * 1 = Probably not sensible!
	 * 2 = Smoothing in slices only
	 * 3 = Smoothing by volume
	 */
    int m_spatial_dims; // 0 = no spatial norm; 2 = slice only; 3 = volume

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

    /**
     * Shrinkage prior type. 
     *
     * One of m, M, p, P or S. Only one of these can be specified in a given run.
     */
    char m_shrinkage_type;

    /**
	 * Nearest-neighbours of each voxel. Vector size is number of voxels,
	 * each entry is vector of indices of nearest neighbours, starting at 1.
	 *
	 * FIXME Sparse matrix would be better?
	 */
    std::vector<std::vector<int> > m_neighbours;

    /**
	 * Next-nearest-neighbours of each voxel. Vector size is number of voxels,
	 * each entry is vector of indices of second nearest neighbours, starting at 1.
	 *
	 * FIXME Sparse matrix would be better?
	 */
    std::vector<std::vector<int> > m_neighbours2;

    // For the new (Sahani-based) smoothing method:
    CovarianceCache m_covar;

    // StS matrix used for S and Z spatial priors
    NEWMAT::SymmetricMatrix m_sts;

    /**
	 * How to measure distances between voxels.
	 *
	 * dist1 = Euclidian distance
	 * dist2 = squared Euclidian distance
	 * mdist = Manhattan distance (|dx| + |dy|)
	 */
    std::string m_dist_measure;

    double m_fixed_delta;
    double m_fixed_rho;

    /**
	 * Update spatial priors on first iteration?
	 */
    bool m_update_first_iter;

    /**
	 * Use evidence optimization
	 */
    bool m_use_evidence;

    double m_always_inital_delta_guess;

    /**
	 * Use full evidence optimization
	 */
    bool m_use_full_evidence;

    /**
	 * Use simultaneous evidence optimization
	 */
    bool m_use_sim_evidence;

    bool m_save_without_prior;
    bool m_save_spatial_priors;
    bool m_locked_linear;

    int m_full_eo_first_param;
    bool m_use_covar_marginals_not_precisions;
    bool m_keep_param_covars;

    int m_new_delta_evaluations;

    bool m_brute_force_delta_search;
};
