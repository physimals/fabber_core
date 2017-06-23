/** 
 * prior.cc
 *
 * Class for a parameter prior
 *
 * Copyright (C) 2007-2017 University of Oxford  
 */

/*  CCOPYRIGHT */

#include "priors.h"

#include "rundata.h"
#include "dist_mvn.h"

#include <newmat.h>
#include <miscmaths/miscmaths.h>

#include <math.h>
#include <vector>
#include <string>
#include <ostream>

using namespace std;
using namespace NEWMAT;
using MISCMATHS::digamma;

std::ostream &operator<<(std::ostream &out, const Prior &prior)
{
    prior.DumpInfo(out);
    return out;
}

DefaultPrior::DefaultPrior(char type, unsigned int idx, string param_name, double mean, double prec)
    : m_param_name(param_name)
    , m_idx(idx)
    , m_type_code(type)
    , m_mean(mean)
    , m_prec(prec)
{
}

void DefaultPrior::DumpInfo(std::ostream &out) const
{
    out << "DefaultPrior: Parameter " << m_idx << " '" << m_param_name << "'" << " mean: "
        << m_mean << " precision: " << m_prec;
}

double DefaultPrior::ApplyToMVN(MVNDist *prior, const PriorContext &ctx)
{
    prior->means(m_idx + 1) = m_mean;

    SymmetricMatrix prec = prior->GetPrecisions();
    prec(m_idx + 1, m_idx + 1) = m_prec;
    prior->SetPrecisions(prec);

    return 0;
}

ImagePrior::ImagePrior(unsigned int idx, string param_name, string filename, double prec, FabberRunData &rundata)
    : DefaultPrior('I', idx, param_name, -1, prec), m_filename(filename)
{
    m_log = rundata.GetLogger();
    m_image = rundata.GetVoxelData(m_filename).AsRow();
}

void ImagePrior::DumpInfo(std::ostream &out) const
{
    out << "ImagePrior: Parameter " << m_idx << " '" << m_param_name << "'" << " filename: "
        << m_filename << " precision: " << m_prec;
}

double ImagePrior::ApplyToMVN(MVNDist *prior, const PriorContext &ctx)
{
    prior->means(m_idx + 1) = m_image(ctx.v);

    SymmetricMatrix prec = prior->GetPrecisions();
    prec(m_idx + 1, m_idx + 1) = m_prec;
    prior->SetPrecisions(prec);

    return 0;
}

void ARDPrior::DumpInfo(std::ostream &out) const
{
    out << "ARDPrior: Parameter " << m_idx << " '" << m_param_name << "'";
}

// Calculate log-gamma from a Taylor expansion; good to one part in 2e-10.
static double gammaln(double x)
{
    ColumnVector series(7);
    series << 2.5066282746310005 << 76.18009172947146 << -86.50532032941677 << 24.01409824083091 << -1.231739572450155
           << 0.1208650973866179e-2 << -0.5395239384953e-5;

    double total = 1.000000000190015;
    for (int i = 2; i <= series.Nrows(); i++)
        total += series(i) / (x + i - 1);

    return log(series(1) * total / x) + (x + 0.5) * log(x + 5.5) - x - 5.5;
}

double ARDPrior::ApplyToMVN(MVNDist *prior, const PriorContext &ctx)
{
    SymmetricMatrix cov = prior->GetCovariance();

    if (ctx.it == 0) {
        // Special case for first iter
        cov(m_idx+1, m_idx+1) = 1e12; //set prior to be initally non-informative
        prior->means(m_idx+1) = 0;
        return 0; //FIXME
    }
    else {
        // Update covariance on subsequent iterations
        // (Chappel et al 2009 Eq D4)
        double post_mean = ctx.fwd_post[ctx.v-1].means(m_idx+1);
        double post_cov = ctx.fwd_post[ctx.v-1].GetCovariance()(m_idx+1, m_idx+1);
        double new_cov = post_mean*post_mean + post_cov;
        cov(m_idx+1, m_idx+1) = new_cov;
        prior->SetCovariance(cov);

        // Calculate the free energy contribution from ARD term
        // (Chappel et al 2009, end of Appendix D)
        double b = 2 / new_cov;
        return -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5 * log(b); 
    }
}

SpatialPrior::SpatialPrior(char type, unsigned int idx, std::string param_name, double mean, double prec, FabberRunData &rundata) 
    : DefaultPrior(type, idx, param_name, mean, prec), m_akmean(1e-8), m_spatial_dims(3), m_spatial_speed(-1) 
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

    //FIXME check valid range
    m_spatial_speed = rundata.GetDoubleDefault("spatial-speed", -1);

    // FIXME still needed?
    m_update_first_iter = rundata.GetBool("update-spatial-prior-on-first-iteration");
}

double SpatialPrior::CalculateAkmean(const PriorContext &ctx)
{
    // The following calculates Tr[Sigmak*S'*S]
    double tmp1 = 0.0;
    double tmp2 = 0.0;
    for (int v = 1; v <= ctx.nvoxels; v++)
    {
        double sigmak = ctx.fwd_post.at(v - 1).GetCovariance()(m_idx+1, m_idx+1);
        int nn = ctx.neighbours.at(v - 1).size();
        if (m_type_code == 'm') // useMRF)
            tmp1 += sigmak * m_spatial_dims * 2;
        else if (m_type_code == 'M') // useMRF2)
            tmp1 += sigmak * (nn + 1e-8);
        else if (m_type_code == 'p')
            tmp1 += sigmak * (4 * m_spatial_dims * m_spatial_dims + nn);
        else // P
            tmp1 += sigmak * (nn*nn + nn);

        double wk = ctx.fwd_post.at(v - 1).means(m_idx+1);
        double Swk = 0.0;
        for (vector<int>::const_iterator v2It = ctx.neighbours[v - 1].begin();
                v2It != ctx.neighbours.at(v - 1).end(); ++v2It)
        {
            Swk += wk - ctx.fwd_post.at(*v2It-1).means(m_idx+1);
        }
        if (m_type_code == 'p' || m_type_code == 'm')
            Swk += wk * (m_spatial_dims * 2 - ctx.neighbours.at(v - 1).size());

        if (m_type_code == 'm' || m_type_code == 'M')
            tmp2 += Swk*wk;
        else 
            tmp2 += Swk*Swk;
    }

    LOG << "SpatialPrior::UpdateAkmean " << m_idx << ": tmp1=" << tmp1 << ", tmp2=" << tmp2 << endl;

    double gk = 1 / (0.5 * tmp1 + 0.5 * tmp2 + 0.1); // prior q1 == 10 (1/q1 == 0.1)
    double akmean = gk * (ctx.nvoxels * 0.5 + 1.0); // prior q2 == 1.0
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
        LOG << "SpatialPrior::UpdateAkmean " << m_idx << ": Rate-limiting the increase on akmean: was "
            << akmean << ", now " << akmeanMax << endl;
        akmean = akmeanMax;
    }

    LOG << "SpatialPrior::UpdateAkmean " << m_idx << ": New akmean: " << akmean << endl;
    return akmean;
}

void SpatialPrior::DumpInfo(std::ostream &out) const
{
    out << "SpatialPrior: Parameter " << m_idx << " '" << m_param_name << "'" << " type " << m_type_code << " mean: "
        << m_mean << " precision: " << m_prec;
}

double SpatialPrior::ApplyToMVN(MVNDist *prior, const PriorContext &ctx)
{
    if (ctx.v == 1 && (ctx.it > 0 || m_update_first_iter)) {
        m_akmean = CalculateAkmean(ctx);
    }

    double weight8 = 0;    // weighted +8
    double contrib8 = 0.0;
    for (vector<int>::const_iterator nidIt = ctx.neighbours[ctx.v - 1].begin();
         nidIt != ctx.neighbours[ctx.v - 1].end(); ++nidIt)
    {
        int nid = *nidIt;
        const MVNDist &neighbourPost = ctx.fwd_post[nid - 1];
        contrib8 += 8 * neighbourPost.means(m_idx+1);
        weight8 += 8;
    }

    double weight12 = 0; // weighted -1, may be duplicated
    double contrib12 = 0.0;
    for (vector<int>::const_iterator nidIt = ctx.neighbours2[ctx.v - 1].begin();
         nidIt != ctx.neighbours2[ctx.v - 1].end(); ++nidIt)
    {
        int nid = *nidIt;
        const MVNDist &neighbourPost = ctx.fwd_post[nid - 1];
        contrib12 += -neighbourPost.means(m_idx+1);
        weight12 += -1;
    }

    int nn = ctx.neighbours[ctx.v - 1].size();

    if (m_type_code == 'p')
    {
        assert(nn <= m_spatial_dims * 2);
        weight8 = 8 * 2 * m_spatial_dims;
        weight12 = -1 * (4 * m_spatial_dims * m_spatial_dims - nn);
    }

    
    double spatial_prec = 0;

    if (m_type_code == 'P')
        spatial_prec = m_akmean * (nn*nn + nn);
    else if (m_type_code == 'm')
        spatial_prec = m_akmean * m_spatial_dims * 2;
    else if (m_type_code == 'M')
        spatial_prec = m_akmean * (nn + 1e-8);
    else if (m_type_code == 'p')
        spatial_prec = m_akmean * (4 * m_spatial_dims * m_spatial_dims + nn);
    else 
        assert(false);

    SymmetricMatrix precs = prior->GetPrecisions();
    if (m_type_code == 'p' || m_type_code == 'm')
    {
        //	Penny-style DirichletBC priors -- ignoring initialFwdPrior completely!
        precs(m_idx+1, m_idx+1) = spatial_prec;
    }
    else
    {
        precs(m_idx+1, m_idx+1) = m_prec + spatial_prec;
    }
    prior->SetPrecisions(precs);

    double mTmp;
    if (m_type_code == 'm') {
        // Dirichlet BCs on MRF
        mTmp = contrib8 / (8 * m_spatial_dims * 2); 
    }               
    else if (m_type_code == 'M') {
        mTmp = contrib8 / (8 * (nn + 1e-8));
    }
    else if (weight8 != 0) {
        mTmp = (contrib8 + contrib12) / (weight8 + weight12);
    }
    else
        mTmp = 0;

    if (m_type_code == 'm' || m_type_code == 'M')
        prior->means(m_idx+1) = prior->GetCovariance()(m_idx+1, m_idx+1) * spatial_prec * mTmp; // = mTmp for p or m
    else {
        // equivalent, when non-spatial priors are very weak: m_fwd_prior[v-1].means = mTmp;
        prior->means(m_idx+1) = prior->GetCovariance()(m_idx+1, m_idx+1) * (spatial_prec * mTmp + m_prec * m_mean);
    }

    return 0;
}

PriorFactory::PriorFactory(const FwdModel &model, FabberRunData &rundata) 
    : m_model(model), m_rundata(rundata)
{
    m_log = rundata.GetLogger();
    m_model.NameParams(m_param_names);
}

vector<Prior *> PriorFactory::CreatePriors()
{
    // Hacky way to get model default priors
    size_t nparams = m_param_names.size();
    MVNDist model_prior(nparams, m_log), model_post(nparams, m_log);
    m_model.HardcodedInitialDists(model_prior, model_post);
    SymmetricMatrix precs = model_prior.GetPrecisions();

    vector<Prior *> priors;
    for (size_t i = 0; i < nparams; i++)
    {
        Prior *prior = PriorFactory::CreatePrior(i, model_prior.means(i+1), precs(i+1, i+1));
        LOG << "PriorFactory::CreatePriors " << *prior << endl;
        priors.push_back(prior);
    }

    return priors;
}

Prior *PriorFactory::CreatePrior(unsigned int idx, double mean, double prec)
{
    char type = m_rundata.GetStringDefault("default-prior-type", "-")[0];
    string img_filename = "";
    
    // Complexity below is due to there being two ways of specifying
    // priors. One is using the param-spatial-priors option which is
    // a sequence of chars in model parameter order, one for each
    // parameter. A + character means 'use the previous value for all
    // remaining parameters'. An 'I' means an image prior and
    // the filename is specified separately using an image-prior<n> option
    string types = GetTypesString();
    char current_type = '-';
    for (size_t i = 0; i < types.size(); i++)
    {
        if (i == idx)
        {
            if (types[i] == '+')
                type = current_type;
            else
                type = types[i];
            break;
        }
        else if (types[i] != '+')
        {
            current_type = types[i];
        }
    }
    // Record the data key (filename) for an image prior. Note that the index is
    // conceptually different from the PSP_byname_image method use below - here
    // it is the parameter index in the model's list (starting at 1), below it depends on
    // the order in which the names are given in the options. Also an optional precision
    // may be set, -1 means unset (since precisions must be > 0)
    if (type == 'I') img_filename = "image-prior" + stringify(idx + 1);

    // Here is the second way of specifying priors which will override the above.
    // PSP_byname<n>=parameter name
    // PSP_byname<n>_type=character indicating type
    // PSP_byname<n>_image=image prior filename
    // PSP_byname<n>_prec=image prior precision
    int current = 0;
    while (true)
    {
        current++;
        string name = m_rundata.GetStringDefault("PSP_byname" + stringify(current), "stop!");
        if (name == "stop!")
            break;

        if (name == m_param_names[idx])
        {
            type = convertTo<char>(m_rundata.GetStringDefault("PSP_byname" + stringify(current) + "_type", stringify(type)));
            if (type == 'I') img_filename = "PSP_byname" + stringify(current) + "_image";

            // Can override mean and precision on command line
            mean = m_rundata.GetDoubleDefault("PSP_byname" + stringify(current) + "_mean", mean);
            prec = m_rundata.GetDoubleDefault("PSP_byname" + stringify(current) + "_prec", prec);
        }
    }

    switch (type) {case 'N':
        case '-':
            return new DefaultPrior(type, idx, m_param_names[idx], mean, prec);
        case 'I':
            return new ImagePrior(idx, m_param_names[idx], img_filename, prec, m_rundata);
        case 'M':
        case 'm':
        case 'P':
        case 'p':
            return new SpatialPrior(type, idx, m_param_names[idx], mean, prec, m_rundata);
        case 'A':
            return new ARDPrior(idx, m_param_names[idx]);
        default:
            throw InvalidOptionValue("DefaultPrior type", stringify(type), "Supported types: NMmPpAI");
    }
}

string PriorFactory::GetTypesString()
{
    size_t num_params = m_param_names.size();
    string priors_str = m_rundata.GetStringDefault("param-spatial-priors", "");
    
    // Yuk
    size_t n_str_params = 0;
    for (size_t i=0; i<priors_str.size(); i++) 
    {
        if (priors_str[i] != '+') n_str_params++;
    }

    if (n_str_params > num_params)
    {
        throw InvalidOptionValue("param-spatial-priors", priors_str, "Too many parameters");
    }

    if (priors_str.size() < num_params)
    {
        // Expand '+' char, if present, to give correct number of parameters
        // If not, don't worry will just use default prior type for the missing
        int deficit = num_params - priors_str.size();
        size_t plus_pos = priors_str.find("+");
        if (plus_pos != std::string::npos)
        {
            priors_str.insert(plus_pos, deficit, '+');
        }
    }

    return priors_str;
}