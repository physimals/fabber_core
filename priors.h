/** 
 * prior.h
 *
 * Classes for parameter priors
 *
 * Copyright (C) 2007-2017 University of Oxford  
 */

#include "rundata.h"
#include "dist_mvn.h"
#include "fwdmodel.h"

#include <newmat.h>

#include <vector>
#include <string>
#include <ostream>

// Temp structure for info required by ApplyToMVN.
//
// Should go into inference method as a kind of generic run context.
struct PriorContext 
{
    PriorContext(int nvoxels, std::vector<MVNDist> &fwd_post, std::vector<std::vector<int> > &neighbours, std::vector<std::vector<int> > &neighbours2)
    : it(0), v(1), nvoxels(nvoxels), fwd_post(fwd_post), neighbours(neighbours), neighbours2(neighbours2) {}
    int it;
    int v;
    int nvoxels;
    std::vector<MVNDist> &fwd_post;
    std::vector<std::vector<int> > &neighbours;
    std::vector<std::vector<int> > &neighbours2;
};

/**
 * Abstract interface for a parameter prior
 */
class Prior : public Loggable
{
public:
    virtual ~Prior() {}

    /** Dump info to output stream */
    virtual void DumpInfo(std::ostream &out) const = 0;

    /** 
     * Apply prior information to an MVN 
     *
     * Returns any additional free energy contribution (e.g. for ARD priors)
     */
    virtual double ApplyToMVN(MVNDist *prior, const PriorContext &ctx) = 0;
};

/**
 * Prior which has a mean and precision
 */
class DefaultPrior : public Prior
{
public:
    DefaultPrior();
    DefaultPrior(char type, unsigned int idx, std::string param_name, double mean=-1, double prec=-1);
    virtual ~DefaultPrior() {}

    /** Parameter name this prior applies to */
    std::string m_param_name;

    /** Parameter index number */
    unsigned int m_idx;

    /** DefaultPrior type code */
    char m_type_code;

    /** DefaultPrior mean */
    double m_mean;

    /** DefaultPrior precision */
    double m_prec;

    virtual void DumpInfo(std::ostream &out) const;
    virtual double ApplyToMVN(MVNDist *prior, const PriorContext &ctx);
};

/**
 * Prior which has takes its mean from a per-voxel image, with a constant precision
 */
class ImagePrior : public DefaultPrior
{
public:
    ImagePrior(unsigned int idx, std::string param_name, std::string filename, double prec, FabberRunData &rundata);

    virtual void DumpInfo(std::ostream &out) const;
    virtual double ApplyToMVN(MVNDist *prior, const PriorContext &ctx);
protected:
    /** Filename containing image data if required */
    std::string m_filename;

    /** Image data if required */
    NEWMAT::RowVector m_image;
};

/**
 * ARD prior
 */
class ARDPrior : public DefaultPrior
{
public:
    ARDPrior(unsigned int idx, std::string param_name) : DefaultPrior('A', idx, param_name) {}

    virtual void DumpInfo(std::ostream &out) const;
    virtual double ApplyToMVN(MVNDist *prior, const PriorContext &ctx);
};

/**
 * Prior which uses spatial information to inform the prior
 *
 * Currently this merges all different spatial priors into one class. Ideally
 * find a way to split each one into a separate subclass.
 */
class SpatialPrior : public DefaultPrior
{
public:
    SpatialPrior(char type, unsigned int idx, std::string param_name, double mean, double prec, FabberRunData &rundata);

    virtual void DumpInfo(std::ostream &out) const;
    virtual double ApplyToMVN(MVNDist *prior, const PriorContext &ctx);
protected:
    double CalculateAkmean(const PriorContext &ctx);
    double m_akmean;
    int m_spatial_dims;
    double m_spatial_speed;
    bool m_update_first_iter;
};

/**
 * Creates instances of Prior depending on the input options
 */
class PriorFactory : public Loggable
{
public:
    PriorFactory(const FwdModel &model, FabberRunData &rundata);

    /** Create priors for all model parameters */
    std::vector<Prior *> CreatePriors();

private:
    const FwdModel &m_model;
    FabberRunData &m_rundata;
    std::vector<std::string> m_param_names;

    /** Create a prior for the parameter at the specified index */
    Prior *CreatePrior(unsigned int idx, double mean, double prec);

    std::string GetTypesString();
};

std::ostream &operator<<(std::ostream &out, const Prior &value);
