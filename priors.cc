/**
 * prior.cc
 *
 * Class for a parameter prior
 *
 * Copyright (C) 2007-2017 University of Oxford
 */

/*  CCOPYRIGHT */

#include "priors.h"

#include "dist_mvn.h"
#include "rundata.h"
#include "tools.h"

#include <miscmaths/miscmaths.h>
#include <newmat.h>

#include <math.h>
#include <ostream>
#include <string>
#include <vector>

using namespace std;
using namespace NEWMAT;
using MISCMATHS::digamma;

std::ostream &operator<<(std::ostream &out, const Prior &prior)
{
    prior.DumpInfo(out);
    return out;
}

string Prior::ExpandPriorTypesString(string priors_str, unsigned int num_params)
{
    // Find out how many prior types are in the string, and what the + character
    // should be interpreted as
    unsigned int n_str_params = 0;
    char repeat_type = '-';
    bool plus_found = false;
    for (size_t i = 0; i < priors_str.size(); i++)
    {
        if (priors_str[i] != '+')
        {
            if (!plus_found)
                repeat_type = priors_str[i];
            n_str_params++;
        }
        else if (plus_found)
        {
            throw InvalidOptionValue(
                "param-spatial-priors", priors_str, "Only one + character allowed");
        }
        else
        {
            plus_found = true;
        }
    }

    if (n_str_params > num_params)
    {
        throw InvalidOptionValue("param-spatial-priors", priors_str, "Too many parameters");
    }
    else if (n_str_params < num_params)
    {
        // Expand '+' char, if present, to give correct number of parameters
        // If there is no +, append with '-', meaning 'model default'
        int deficit = num_params - n_str_params;
        size_t plus_pos = priors_str.find("+");
        if (plus_pos != std::string::npos)
        {
            priors_str.insert(plus_pos, deficit - 1, '+');
        }
        else
        {
            priors_str.insert(priors_str.end(), deficit, '-');
        }
    }
    else
    {
        // We already have enough types for all the parameters so erase any
        // pointless + char
        priors_str.erase(std::remove(priors_str.begin(), priors_str.end(), '+'), priors_str.end());
    }

    // Finally, replace all + chars with identified repeat type
    std::replace(priors_str.begin(), priors_str.end(), '+', repeat_type);
    assert(priors_str.size() == num_params);

    return priors_str;
}

DefaultPrior::DefaultPrior(const Parameter &p)
    : m_param_name(p.name)
    , m_idx(p.idx)
    , m_type_code(p.prior_type)
    , m_params(p.prior)
{
}

void DefaultPrior::DumpInfo(std::ostream &out) const
{
    out << "DefaultPrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " mean: " << m_params.mean() << " precision: " << m_params.prec();
}

double DefaultPrior::ApplyToMVN(MVNDist *prior, const ThreadContext &ctx)
{
    prior->means(m_idx + 1) = m_params.mean();

    SymmetricMatrix prec = prior->GetPrecisions();
    prec(m_idx + 1, m_idx + 1) = m_params.prec();
    prior->SetPrecisions(prec);

    return 0;
}

ImagePrior::ImagePrior(const Parameter &p, FabberRunData &rundata)
    : DefaultPrior(p)
{
    m_log = rundata.GetLogger();
    m_filename = p.options.find("image")->second;
    m_image = rundata.GetVoxelData(m_filename).AsRow();
}

void ImagePrior::DumpInfo(std::ostream &out) const
{
    out << "ImagePrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " filename: " << m_filename << " precision: " << m_params.prec();
}

double ImagePrior::ApplyToMVN(MVNDist *prior, const ThreadContext &ctx)
{
    prior->means(m_idx + 1) = m_image(ctx.v);

    SymmetricMatrix prec = prior->GetPrecisions();
    prec(m_idx + 1, m_idx + 1) = m_params.prec();
    prior->SetPrecisions(prec);

    return 0;
}

void ARDPrior::DumpInfo(std::ostream &out) const
{
    out << "ARDPrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " initial mean: " << m_params.mean() << " initial precision: " << m_params.prec();
}

double ARDPrior::ApplyToMVN(MVNDist *prior, const ThreadContext &ctx)
{
    SymmetricMatrix cov = prior->GetCovariance();
    double post_mean = ctx.fwd_post[ctx.v - 1].means(m_idx + 1);
    double post_cov = ctx.fwd_post[ctx.v - 1].GetCovariance()(m_idx + 1, m_idx + 1);
    // (Chappel et al 2009 Eq D4)
    double new_cov = post_mean * post_mean + post_cov;

    if (ctx.it == 0)
    {
        // Special case for first iter
        // Initially set prior to model default. The other alternative is to set
        // it to be initially non-informative, however this can still be achieved
        // by specifying the model prior mean/precision by options.
        cov(m_idx + 1, m_idx + 1) = m_params.var();
        prior->means(m_idx + 1) = m_params.mean();
        // LOG << "first iter ARD: " << m_params.var() << ", " << m_params.mean() << endl;
    }
    else
    {
        // LOG << "post: " << post_cov << ", " << post_mean << endl;
        // Update covariance on subsequent iterations
        cov(m_idx + 1, m_idx + 1) = new_cov;
        // LOG << "subs iter ARD: " << new_cov << ", " << prior->means(m_idx+1) << endl;
    }
    prior->SetCovariance(cov);

    // Calculate the free energy contribution from ARD term
    // (Chappel et al 2009, end of Appendix D)
    double b = 2 / new_cov;
    return -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5 * log(b);
}

SpatialPrior::SpatialPrior(const Parameter &p, FabberRunData &rundata)
    : DefaultPrior(p)
    , m_akmean(1e-8)
    , m_spatial_dims(3)
    , m_spatial_speed(-1)
{
    m_log = rundata.GetLogger();
    m_spatial_dims = rundata.GetIntDefault("spatial-dims", 3);
    if (m_spatial_dims < 0 || m_spatial_dims > 3)
    {
        throw InvalidOptionValue("spatial-dims", stringify(m_spatial_dims), "Must be 0, 1, 2 or 3");
    }
    else if (m_spatial_dims == 1)
    {
        WARN_ONCE("spatial-dims=1 is very weird... hope you're just testing!");
    }
    else if (m_spatial_dims == 2)
    {
        WARN_ONCE("spatial-dims=2 doesn't decompose into slices");
    }

    // FIXME check valid range
    m_spatial_speed = rundata.GetDoubleDefault("spatial-speed", -1);

    // FIXME still needed?
    m_update_first_iter = rundata.GetBool("update-spatial-prior-on-first-iteration");
}

double SpatialPrior::CalculateAkmean(const SpatialVbThreadContext &sctx)
{
    // The following calculates Tr[Sigmak*S'*S]
    double tmp1 = 0.0;
    double tmp2 = 0.0;
    for (int v = 1; v <= sctx.nvoxels; v++)
    {
        // Ignore voxels where numerical issues have occurred
        if (std::find(sctx.ignore_voxels.begin(), sctx.ignore_voxels.end(), v)
            != sctx.ignore_voxels.end()) continue;

        double sigmak = sctx.fwd_post.at(v - 1).GetCovariance()(m_idx + 1, m_idx + 1);
        int nn = sctx.m_neighbours->at(v - 1).size();
        if (m_type_code == PRIOR_SPATIAL_m) // useMRF)
            tmp1 += sigmak * m_spatial_dims * 2;
        else if (m_type_code == PRIOR_SPATIAL_M) // useMRF2)
            tmp1 += sigmak * (nn + 1e-8);
        else if (m_type_code == PRIOR_SPATIAL_p)
            tmp1 += sigmak * (4 * m_spatial_dims * m_spatial_dims + nn);
        else // P
            tmp1 += sigmak * (nn * nn + nn);

        double wk = sctx.fwd_post.at(v - 1).means(m_idx + 1);
        double Swk = 0.0;
        for (vector<int>::const_iterator v2It = (*sctx.m_neighbours2)[v - 1].begin();
             v2It != sctx.m_neighbours->at(v - 1).end(); ++v2It)
        {
            Swk += wk - sctx.fwd_post.at(*v2It - 1).means(m_idx + 1);
        }
        if (m_type_code == PRIOR_SPATIAL_p || m_type_code == PRIOR_SPATIAL_m)
            Swk += wk * (m_spatial_dims * 2 - sctx.m_neighbours->at(v - 1).size());

        if (m_type_code == PRIOR_SPATIAL_m || m_type_code == PRIOR_SPATIAL_M)
            tmp2 += Swk * wk;
        else
            tmp2 += Swk * Swk;
    }

    LOG << "SpatialPrior::UpdateAkmean " << m_idx << ": tmp1=" << tmp1 << ", tmp2=" << tmp2 << endl;

    double gk = 1 / (0.5 * tmp1 + 0.5 * tmp2 + 0.1); // prior q1 == 10 (1/q1 == 0.1)
    double akmean = gk * (sctx.nvoxels * 0.5 + 1.0);  // prior q2 == 1.0
    double akmeanMax = akmean * m_spatial_speed;
    if (akmean < 1e-50)
    {
        LOG << "SpatialPrior::UpdateAkmean " << m_idx << ": was " << akmean << endl;
        WARN_ONCE("SpatialPrior::UpdateAkmean akmean value was tiny!");
        akmean = 1e-50; // prevent crashes
    }

    if (akmeanMax < 0.5)
    {
        akmeanMax = 0.5; // totally the wrong scaling.. oh well
    }

    if (m_spatial_speed > 0 && akmean > akmeanMax)
    {
        LOG << "SpatialPrior::UpdateAkmean " << m_idx
            << ": Rate-limiting the increase on akmean: was " << akmean << ", now " << akmeanMax
            << endl;
        akmean = akmeanMax;
    }

    LOG << "SpatialPrior::UpdateAkmean " << m_idx << ": New akmean: " << akmean << endl;
    return akmean;
}

void SpatialPrior::DumpInfo(std::ostream &out) const
{
    out << "SpatialPrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " type " << m_type_code << " mean: " << m_params.mean()
        << " precision: " << m_params.prec();
}

double SpatialPrior::ApplyToMVN(MVNDist *prior, const ThreadContext &ctx)
{
    const SpatialVbThreadContext& sctx = dynamic_cast<const SpatialVbThreadContext&>(ctx);

    if (sctx.v == 1 && (sctx.it > 0 || m_update_first_iter))
    {
        m_akmean = CalculateAkmean(sctx);
    }

    double weight8 = 0; // weighted +8
    double contrib8 = 0.0;
    for (vector<int>::const_iterator nidIt = (*sctx.m_neighbours2)[sctx.v - 1].begin();
         nidIt != (*sctx.m_neighbours2)[sctx.v - 1].end(); ++nidIt)
    {
        int nid = *nidIt;
        const MVNDist &neighbourPost = sctx.fwd_post[nid - 1];
        contrib8 += 8 * neighbourPost.means(m_idx + 1);
        weight8 += 8;
    }

    double weight12 = 0; // weighted -1, may be duplicated
    double contrib12 = 0.0;
    for (vector<int>::const_iterator nidIt = (*sctx.m_neighbours2)[sctx.v - 1].begin();
         nidIt != (*sctx.m_neighbours2)[sctx.v - 1].end(); ++nidIt)
    {
        int nid = *nidIt;
        const MVNDist &neighbourPost = sctx.fwd_post[nid - 1];
        contrib12 += -neighbourPost.means(m_idx + 1);
        weight12 += -1;
    }

    int nn = (*sctx.m_neighbours2)[sctx.v - 1].size();

    if (m_type_code == PRIOR_SPATIAL_p)
    {
        assert(nn <= m_spatial_dims * 2);
        weight8 = 8 * 2 * m_spatial_dims;
        weight12 = -1 * (4 * m_spatial_dims * m_spatial_dims - nn);
    }

    double spatial_prec = 0;

    if (m_type_code == PRIOR_SPATIAL_P)
        spatial_prec = m_akmean * (nn * nn + nn);
    else if (m_type_code == PRIOR_SPATIAL_m)
        spatial_prec = m_akmean * m_spatial_dims * 2;
    else if (m_type_code == PRIOR_SPATIAL_M)
        spatial_prec = m_akmean * (nn + 1e-8);
    else if (m_type_code == PRIOR_SPATIAL_p)
        spatial_prec = m_akmean * (4 * m_spatial_dims * m_spatial_dims + nn);
    else
        assert(false);

    // Set the prior precision for this parameter
    SymmetricMatrix precs = prior->GetPrecisions();
    if (m_type_code == PRIOR_SPATIAL_p || m_type_code == PRIOR_SPATIAL_m)
    {
        //	Penny-style DirichletBC priors -- ignoring initialFwdPrior completely!
        precs(m_idx + 1, m_idx + 1) = spatial_prec;
    }
    else
    {
        precs(m_idx + 1, m_idx + 1) = m_params.prec() + spatial_prec;
    }
    prior->SetPrecisions(precs);

    // Set the prior mean for this parameter
    // Note that we multiply by reciprocals rather than dividing. This is
    // to maximise numerical compatibility with NEWMAT which presumably
    // does it as an optimization when dividing a whole matrix by a constant
    double mTmp;
    if (m_type_code == PRIOR_SPATIAL_m)
    {
        // Dirichlet BCs on MRF
        double rec = 1 / (8 * m_spatial_dims * 2);
        mTmp = contrib8 * rec;
    }
    else if (m_type_code == PRIOR_SPATIAL_M)
    {
        double rec = 1 / (8 * (double(nn) + 1e-8));
        mTmp = contrib8 * rec;
    }
    else if (weight8 != 0)
    {
        double rec = 1 / (weight8 + weight12);
        mTmp = (contrib8 + contrib12) * rec;
    }
    else
        mTmp = 0;

    // LOG << "SpatialPrior:: " << prior->GetCovariance()(m_idx+1, m_idx+1) << ", " << spatial_prec
    // << ", " << contrib8 << ", " << den << ", " << mTmp << " : " << t1 << endl;

    if (m_type_code == PRIOR_SPATIAL_m || m_type_code == PRIOR_SPATIAL_M)
        prior->means(m_idx + 1) = prior->GetCovariance()(m_idx + 1, m_idx + 1) * spatial_prec
            * mTmp; // = mTmp for p or m
    else
    {
        // equivalent, when non-spatial priors are very weak: m_fwd_prior[v-1].means = mTmp;
        prior->means(m_idx + 1) = prior->GetCovariance()(m_idx + 1, m_idx + 1)
            * (spatial_prec * mTmp + m_params.prec() * m_params.mean());
    }

    return 0;
}

std::vector<Prior *> PriorFactory::CreatePriors(const std::vector<Parameter> &params)
{
    vector<Prior *> priors;
    for (size_t i = 0; i < params.size(); i++)
    {
        Prior *prior = PriorFactory::CreatePrior(params[i]);
        LOG << "PriorFactory::CreatePriors " << *prior << endl;
        priors.push_back(prior);
    }

    return priors;
}

PriorFactory::PriorFactory(FabberRunData &rundata)
    : Loggable(rundata.GetLogger())
    , m_rundata(rundata)
{
}

Prior *PriorFactory::CreatePrior(Parameter p)
{
    switch (p.prior_type)
    {
    case PRIOR_NORMAL:
    case PRIOR_DEFAULT:
        return new DefaultPrior(p);
    case PRIOR_IMAGE:
        return new ImagePrior(p, m_rundata);
    case PRIOR_SPATIAL_M:
    case PRIOR_SPATIAL_m:
    case PRIOR_SPATIAL_P:
    case PRIOR_SPATIAL_p:
        return new SpatialPrior(p, m_rundata);
    case PRIOR_ARD:
        return new ARDPrior(p, m_rundata);
    default:
        throw InvalidOptionValue("Prior type", stringify(p.prior_type), "Supported types: NMmPpAI");
    }
}
