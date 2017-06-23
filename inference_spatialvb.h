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
        , m_update_first_iter(false)
        , m_locked_linear(false)
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
    void SetFwdPriorShrinkageType(int voxel, const NEWMAT::DiagonalMatrix &akmean);
    double SetFwdPrior(int voxel, bool isFirstIteration);

    /**
	 * Calculate first and second nearest neighbours of each voxel
	*/
    void CalcNeighbours(const NEWMAT::Matrix &voxelCoords);

    /**
     * Ignore this voxel in future updates.
     *
     * No calculation of priors or posteriors will occur for this voxel
     * and it will be removed from the lists of neighbours for other voxels.
     * The effect should be as if it were masked
     */
    void IgnoreVoxel(int v);

    // Per-voxel prior and posterior distributions. For Spatial VB we need to
    // keep these around during iteration as the influence the calculations on
    // neighbouring voxels. In particular the priors change as they reflect
    // the assumption on spatial smoothness
    //
    // NoiseParams is polymorphic type, so need to use pointers
    std::vector<NoiseParams *> m_noise_post;
    std::vector<NoiseParams *> m_noise_prior;
    std::vector<MVNDist> m_fwd_prior;
    std::vector<MVNDist> m_fwd_post;
    std::vector<LinearizedFwdModel> m_lin_model;

    /**
     * Voxels to ignore, indexed from 1 as per NEWMAT
     */
    std::vector<int> m_ignore_voxels;
    
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
	 * Maximum precision increase per iteration
     * 
     * If >0 but <1, restrict akmeans by this factor
     * If >1 or <0, no restriction
	 */
    double m_spatial_speed;

    /**
	 * Type of spatial prior to use for each parameter. Should be one
	 * character per parameter, however if string ends with + then
	 * the last character is repeated for remaining parameters
	 */
    std::string m_prior_types_str;

    /**
     * Shrinkage prior type. 
     *
     * One of m, M, p, P. Only one of these can be specified in a given run.
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

    /**
	 * Update spatial priors on first iteration?
	 */
    bool m_update_first_iter;

    /**
     * Reduce this to a linear problem, using the given
     * voxelwise linearizations (probably loaded from an MVN)
     */
    std::string m_locked_linear_file;

    bool m_locked_linear;
};
