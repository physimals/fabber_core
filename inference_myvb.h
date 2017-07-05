/*  inference_spatialvb.h - implementation of VB with spatial priors

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference.h"
#include "convergence.h"
#include "run_context.h"

#include <string>
#include <vector>

class Vb : public InferenceTechnique
{
public:
    static InferenceTechnique *NewInstance();

    virtual void GetOptions(vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;
    virtual string GetVersion() const;

    Vb()
        : m_nvoxels(0)
        , m_noise_params(0)
        , m_needF(false)
        , m_printF(false)
        , continueFwdOnly(false)
        , m_origdata(NULL)
        , m_coords(NULL)
        , m_suppdata(NULL)
        , m_num_mcsteps(0)
        , m_conv(NULL)
        , m_spatial_dims(-1)
        , m_locked_linear(false)
    {
    }

    virtual void Initialize(FwdModel *fwd_model, FabberRunData &args);
    virtual void DoCalculations(FabberRunData &data);
    //    virtual ~Vb();

    void SaveResults(FabberRunData &rundata) const;
protected:

    std::string GetPriorTypesString(FabberRunData &rundata);
    void InitializeNoiseFromParam(FabberRunData &args, NoiseParams *dist, string param_key);
    void PassModelData(int voxel);

    virtual void DoCalculationsVoxelwise(FabberRunData &data);
    virtual void DoCalculationsSpatial(FabberRunData &data);

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

    // These are used for resuming a previous calculation
    /** If not empty, load initial MVN from this file */
    std::string m_continueFromFile;

    std::string paramFilename;

    /** If set, initial MVN only has fwd model information, not noise */
    bool continueFwdOnly;

    /** Free energy for each voxel */
    std::vector<double> resultFs;

    const NEWMAT::Matrix *m_origdata;
    const NEWMAT::Matrix *m_coords;
    const NEWMAT::Matrix *m_suppdata;

    /** Number of motion correction steps to run */
    int m_num_mcsteps;

    /**
     * Convergence detector in use
     */
    ConvergenceDetector *m_conv;
    
    RunContext *m_ctx;
    
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
