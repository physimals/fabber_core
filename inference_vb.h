/*  inference_spatialvb.h - implementation of VB with spatial priors

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */
#pragma once

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
    virtual void SaveResults(FabberRunData &rundata) const;
protected:
    /**
     * Determine whether we need spatial VB mode
     *
     * It is required either because it has been asked for (--method=spatialvb) or
     * if any spatial priors have been specified (types mMpP)
     */
    bool IsSpatial(ThreadContext *ctx, FabberRunData &rundata) const;

    /**
    * Check voxels are listed in order
    *
    * Order must be increasing in z value, or if same
    * increasing in y value, and if y and z are same
    * increasing in x value.
    *
    * This is basically column-major (Fortran) ordering - used as default by NEWIMAGE.
    */
    void CheckCoordMatrixCorrectlyOrdered(const NEWMAT::Matrix &coords);
    void CalcNeighbours(const NEWMAT::Matrix &coords, int spatial_dims, std::vector<std::vector<int> > &neighbours, std::vector<std::vector<int> > &neighbours2);

    int m_num_mcsteps;
};

class VbThreadContext : public ThreadContext
{
public:
    VbThreadContext(FabberRunData &rundata, int worker_id=0, int n_workers=1, int start_vox=1, int num_vox=-1);

protected:
    /**
     * Do calculations loop in voxelwise mode (i.e. all iterations for
     * one voxel, then all iterations for the next voxel, etc)
     */
    void Run();

    bool m_printF, m_needF, m_saveF;

    /**
     * Calculate free energy if required, and display if required
     */
    double CalculateF(int v, std::string label, double Fprior);

    /**
     * Output detailed debugging information for a voxel
     */
    void DebugVoxel(int v, const string &where);
};

class SpatialVbThreadContext : public VbThreadContext
{
public:
    SpatialVbThreadContext(FabberRunData &rundata, std::vector<std::vector<int> > &neighbours, std::vector<std::vector<int> > &neighbours2, int worker_id=0, int n_workers=1, int start_vox=1, int num_vox=-1)
      : VbThreadContext(rundata, worker_id, n_workers, start_vox, num_vox)
      , m_neighbours(&neighbours)
      , m_neighbours2(&neighbours2)
    {}

    const std::vector<std::vector<int> > *m_neighbours;
    const std::vector<std::vector<int> > *m_neighbours2;

protected:
    void Run();

};
