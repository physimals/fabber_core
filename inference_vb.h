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
        : m_num_mcsteps(0)
    {
    }

    virtual void GetOptions(vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;
    virtual string GetVersion() const;

    virtual void Initialize(FabberRunData &args);
    virtual void DoCalculations(FabberRunData &data);
    /**
     * Do calculations loop in spatial mode (i.e. one iteration of all
     * voxels, then next iteration of all voxels, etc)
     */
    virtual void DoCalculationsSpatial(ThreadContext *ctx, FabberRunData &data);
    virtual void SaveResults(FabberRunData &rundata) const;
protected:
    /**
     * Calculate free energy if required, and display if required
     */
    double CalculateF(ThreadContext *ctx, int v, std::string label, double Fprior);

    /**
     * Output detailed debugging information for a voxel
     */
    void DebugVoxel(ThreadContext *ctx, int v, const string &where);

    /**
     * Determine whether we need spatial VB mode
     *
     * It is required either because it has been asked for (--method=spatialvb) or
     * if any spatial priors have been specified (types mMpP)
     */
    bool IsSpatial(ThreadContext *ctx, FabberRunData &rundata) const;

    int m_num_mcsteps;
};

class VbThreadContext : public ThreadContext
{
public:
    VbThreadContext(FabberRunData &rundata, int worker_id=0, int n_workers=1, int start_vox=1, int num_vox=-1)
      : ThreadContext(rundata, worker_id, n_workers, start_vox, num_vox)
    {}
    
    static void *Run(void *obj)
    {
        VbThreadContext *vbc  = ((VbThreadContext *)obj);
        vbc->DoCalculationsVoxelwise();
        return NULL;
    }

    /**
     * Do calculations loop in voxelwise mode (i.e. all iterations for
     * one voxel, then all iterations for the next voxel, etc)
     */
    void DoCalculationsVoxelwise();

    /**
     * Calculate free energy if required, and display if required
     */
    double CalculateF(int v, std::string label, double Fprior);

    /**
     * Output detailed debugging information for a voxel
     */
    void DebugVoxel(int v, const string &where);
};
