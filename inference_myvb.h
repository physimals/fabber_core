/*  inference_spatialvb.h - implementation of VB with spatial priors

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_vb.h"

#include "covariance_cache.h"

class Vb : public VariationalBayesInferenceTechnique
{
public:
    static InferenceTechnique *NewInstance();

    virtual void GetOptions(vector<OptionSpec> &opts) const;

    Vb()
        : VariationalBayesInferenceTechnique()
        , m_spatial_dims(-1)
        , m_locked_linear(false)
    {
    }

    virtual void Initialize(FwdModel *fwd_model, FabberRunData &args);
    virtual void DoCalculations(FabberRunData &data);
    //    virtual ~Vb();

    void SaveResults(FabberRunData &rundata) const;
protected:

    virtual void DoCalculationsVoxelwise(FabberRunData &data);
    virtual void DoCalculationsSpatial(FabberRunData &data);

    std::string GetPriorTypesString(FabberRunData &rundata);

    double CalculateF(int v, std::string label, double Fprior);

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

    std::vector<NoiseParams *> m_noise_post;
    std::vector<NoiseParams *> m_noise_prior;
    std::vector<MVNDist> m_fwd_prior;
    std::vector<MVNDist> m_fwd_post;
    std::vector<LinearizedFwdModel> m_lin_model;
    std::vector<std::vector<int> > m_neighbours;
    std::vector<std::vector<int> > m_neighbours2;
    
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
	 * Type of spatial prior to use for each parameter. Should be one
	 * character per parameter, however if string ends with + then
	 * the last character is repeated for remaining parameters
	 */
    std::string m_prior_types_str;

    /**
     * Reduce this to a linear problem, using the given
     * voxelwise linearizations (probably loaded from an MVN)
     */
    std::string m_locked_linear_file;

    bool m_locked_linear;
};
