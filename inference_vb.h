/*  inference_spatialvb.h - implementation of VB with spatial priors

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "convergence.h"
#include "inference.h"
#include "run_context.h"

#include <string>
#include <vector>

class Vb : public InferenceTechnique
{
public:
    static InferenceTechnique *NewInstance();

    Vb()
        : m_needF(false)
        , m_printF(false)
        , m_saveF(false)
        , m_num_mcsteps(0)
        , m_spatial_dims(-1)
        , m_locked_linear(false)
    {
    }

    virtual void GetOptions(vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;
    virtual string GetVersion() const;

    virtual void Initialize(FabberRunData &args);
    virtual void DoCalculations(FabberRunData &data);

    virtual void SaveResults(FabberRunData &rundata) const;

protected:
    /**
     * Determine whether we need spatial VB mode
     *
     * It is required either because it has been asked for (--method=spatialvb) or
     * if any spatial priors have been specified (types mMpP)
     */
    bool IsSpatial(FabberRunData &rundata) const;

    /**
     * Do calculations loop in voxelwise mode (i.e. all iterations for
     * one voxel, then all iterations for the next voxel, etc)
     */
    virtual void DoCalculationsVoxelwise(FabberRunData &data);

    /**
     * Do calculations loop in spatial mode (i.e. one iteration of all
     * voxels, then next iteration of all voxels, etc)
     */
    virtual void DoCalculationsSpatial(FabberRunData &data);

    /**
     * Calculate free energy if required, and display if required
     */
    double CalculateF(int v, std::string label, double Fprior);

    /**
     * Output detailed debugging information for a voxel
     */
    void DebugVoxel(int v, const string &where);

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

    /** True if convergence detector requires the free energy */
    bool m_needF;

    /** True if we need to print the free energy at each iteration */
    bool m_printF;

    /** True if we need to to save the final free energy */
    bool m_saveF;

    /** Number of motion correction steps to run */
    int m_num_mcsteps;

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
     * Fix the linearization centres of the linearized forward model.
     *
     * This reduces the inference to a purely linear problem. The fixed
     * centres are generally loaded from an MVN file
     */
    bool m_locked_linear;
};
