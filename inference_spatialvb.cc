/*  inference_spatialvb.cc - implementation of VB with spatial priors

 Adrian Groves and Matthew Webster, FMRIB Image Analysis Group

 Copyright (C) 2007-2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_spatialvb.h"

#include "convergence.h"
#include "easylog.h"
#include "tools.h"

#include <miscmaths/miscmaths.h>
#include <newmatio.h>

#include <math.h>

using MISCMATHS::sign;

#define NOCACHE 1

static OptionSpec OPTIONS[] = {
    { "spatial-dims", OPT_INT, "Number of spatial dimensions", OPT_NONREQ, "3" },
    { "spatial-speed", OPT_STR, "Number of spatial dimensions", OPT_NONREQ,
        "-1" },
    { "distance-measure", OPT_STR, "", OPT_NONREQ, "dist1" },
    { "param-spatial-priors", OPT_STR,
        "Type of spatial priors for each parameter, as a sequence of characters. "
        "S=spatial, N=nonspatial, D=Gaussian-process-based ",
        OPT_NONREQ, "S+" },
    { "fixed-delta", OPT_STR, "", OPT_NONREQ, "-1" },
    { "fixed-rho", OPT_STR, "", OPT_NONREQ, "0" },
    { "new-delta-iterations", OPT_INT, "", OPT_NONREQ, "10" },
    { "always-initial-delta-guess", OPT_STR, "", OPT_NONREQ, "-1" },
    { "update-spatial-prior-on-first-iteration", OPT_BOOL, "", OPT_NONREQ, "" },
    { "use-simultaneous-evidence-optimization", OPT_BOOL, "", OPT_NONREQ, "" },
    { "use-full-evidence-optimization", OPT_BOOL, "", OPT_NONREQ, "" },
    { "use-evidence-optimization", OPT_BOOL, "", OPT_NONREQ, "" },
    { "use-covariance-marginals", OPT_BOOL, "", OPT_NONREQ, "" },
    { "keep-interparameter-covariances", OPT_BOOL, "", OPT_NONREQ, "" },
    { "brute-force-delta-search", OPT_BOOL, "", OPT_NONREQ, "" },
    { "no-eo", OPT_BOOL, "", OPT_NONREQ, "" },
    { "slow-eo", OPT_BOOL, "", OPT_NONREQ, "" },
    { "" },
};

void SpatialVariationalBayes::GetOptions(vector<OptionSpec> &opts) const
{
    VariationalBayesInferenceTechnique::GetOptions(opts);

    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

InferenceTechnique *SpatialVariationalBayes::NewInstance()
{
    return new SpatialVariationalBayes();
}

/**
 * How many spatial dimensions have been requested?
 */
static int GetSpatialDims(FabberRunData &args)
{
    int dims = args.GetIntDefault("spatial-dims", 3);

    if (dims < 0 || dims > 3)
    {
        throw InvalidOptionValue("spatial-dims", stringify(dims),
            "Must be 0, 1, 2 or 3");
    }
    else if (dims == 1)
    {
        // WARN_ONCE("spatial-dims=1 is very weird... I hope you're just testing
        // something!");
    }
    else if (dims == 2)
    {
        // WARN_ONCE("spatial-dims=2 doesn't decompose into slices and won't help if
        // you're using the D prior");
    }

    return dims;
}

static string GetPriorTypesStr(vector<PriorType> &priors)
{
    string types;
    for (size_t i = 0; i < priors.size(); i++)
    {
        types += priors[i].m_type;
    }
    return types + "\0";
}

void SpatialVariationalBayes::Initialize(FwdModel *fwd_model,
    FabberRunData &args)
{
    // In spatial VB we want to default to spatial priors so we set this parameter
    // if it's not already set by the user. Note that we are doing this before
    // initializing the VB method so it acts as a default
    args.Set("default-prior-type", "S");
    VariationalBayesInferenceTechnique::Initialize(fwd_model, args);

    m_prior_types_str = GetPriorTypesStr(m_prior_types);
    m_spatial_dims = GetSpatialDims(args);
    m_spatial_speed = args.GetDoubleDefault("spatial-speed", -1);
    assert(m_spatial_speed > 1 || m_spatial_speed == -1);

    m_dist_measure = args.GetStringDefault("distance-measure", "dist1");

    // m_shrinkage_type = GetShrinkageType();

    // Which shrinkage prior to use?  For various reasons we can't
    // yet mix different shinkage priors for different parameters (and I can't
    // see any good reason to implement it).
    m_shrinkage_type = '-';
    for (int k = 0; k < m_num_params; k++)
    {
        char type = m_prior_types[k].m_type;
        switch (type)
        {
        case '-':
        case 'R':
        case 'D':
        case 'N':
        case 'F':
        case 'I':
        case 'A':
            break;

        case 'm':
        case 'M':
        case 'p':
        case 'P':
        case 'S':
            if (type != m_shrinkage_type && m_shrinkage_type != '-')
                throw FabberRunDataError(
                    "Sorry, only one type of shrinkage prior at a time, please!\n");
            m_shrinkage_type = type;
            break;

        default:
            throw InvalidOptionValue("param-spatial-priors", stringify(type),
                "Unrecognized spatial prior type");
        }
    }

    // Some unsupported options:
    m_fixed_delta = args.GetDoubleDefault("fixed-delta", -1);
    m_fixed_rho = args.GetDoubleDefault("fixed-rho", 0);
    m_update_first_iter = args.GetBool("update-spatial-prior-on-first-iteration");
    m_keep_param_covars = args.GetIntDefault("new-delta-iterations", 10);
    assert(m_keep_param_covars > 0);

    // Some deprecated options:
    m_use_sim_evidence = args.GetBool("use-simultaneous-evidence-optimization");
    m_use_full_evidence = m_use_sim_evidence || args.GetBool("use-full-evidence-optimization");
    // WARNING: May need to be set to a sensible value in other circumstances!
    m_full_eo_first_param = m_use_full_evidence ? args.GetIntDefault("first-parameter-for-full-eo", 1)
                                                : -999;
    m_use_evidence = m_use_full_evidence || args.GetBool("use-evidence-optimization");
    m_use_covar_marginals_not_precisions = m_use_full_evidence && args.GetBool("use-covariance-marginals");
    m_keep_param_covars = m_use_full_evidence && args.GetBool("keep-interparameter-covariances");
    m_always_inital_delta_guess = args.GetDoubleDefault("always-initial-delta-guess", -1);
    assert(!(m_update_first_iter && !m_use_evidence)); // currently doesn't work, but fixable
    m_brute_force_delta_search = args.GetBool("brute-force-delta-search");

    m_save_without_prior = m_use_evidence; // or other reasons?
    m_save_spatial_priors = false;
    WARN_ONCE("SpatialVariationalBayes::Initialize Not saving finalSpatialPriors.nii.gz -- too huge!!");

    // Locked linearizations, if requested
    m_locked_linear = (lockedLinearFile != "");

    // Preferred way of using these options
    if (!m_use_full_evidence && !args.GetBool("no-eo") && m_prior_types_str.find_first_of("DR") != string::npos)
    {
        m_use_full_evidence = true;
        m_use_evidence = true;
        // m_keep_param_covars = true; // hacky
        m_use_sim_evidence = args.GetBool("slow-eo");
        if (!m_use_sim_evidence)
            WARN_ONCE("SpatialVariationalBayes::Initialize Defaulting to Full (non-simultaneous) Evidence Optimization");
    }

    if (m_prior_types_str.find("F") != string::npos) // F found
    {
        if (m_fixed_delta < 0)
            throw InvalidOptionValue(
                "param-spatial-priors", "F",
                "Must specify fixed-delta > 0 for this type of spatial prior");
    }
    else
    {
        if (m_fixed_delta == -1)
            m_fixed_delta = 0.5; // Default initial value (in mm!)
    }
}

void SpatialVariationalBayes::SetupPerVoxelDists(FabberRunData &allData)
{
    // Initialized in voxel loop below (from file or default as required)
    m_noise_post.resize(m_nvoxels, NULL);
    m_noise_prior.resize(m_nvoxels, NULL);
    m_fwd_post.resize(m_nvoxels);

    // Initialized in voxel loop if required
    m_fwd_post_no_prior.resize(m_nvoxels, NULL);

    // Re-centred in voxel loop below
    m_lin_model.resize(m_nvoxels, LinearizedFwdModel(m_model));

    // Static initialization for all voxels currently
    m_fwd_prior.resize(m_nvoxels, *initialFwdPrior);

    // Loaded from file if required, otherwise initialized during calculation
    resultMVNs.resize(m_nvoxels, NULL);

    // Initialized during calculation
    resultFs.resize(m_nvoxels, 9999); // 9999 is a garbage default value

    // Empty by default. Resized and initialized below if required
    vector<MVNDist *> lockedLinearDists;
    Matrix lockedLinearCentres;

    const int nNoiseParams = initialNoisePrior->OutputAsMVN().GetSize();

    if (m_locked_linear)
    {
        LOG << "SpatialVariationalBayes::Loading fixed linearization centres from the MVN '"
            << lockedLinearFile << "'\nNOTE: This does not check if the correct "
                                   "number of parameters is present!\n";
        MVNDist::Load(lockedLinearDists, lockedLinearFile, allData, m_log);
        lockedLinearCentres.ReSize(m_num_params, m_nvoxels);
    }

    if (m_continueFromFile != "")
    {
        LOG << "SpatialVbInferenceTechnique::Continuing from file "
            << m_continueFromFile << endl;
        InitMVNFromFile(m_continueFromFile, allData, paramFilename);
    }

    if (m_save_without_prior)
    {
        resultMVNsWithoutPrior.resize(m_nvoxels, NULL);
    }

    for (int v = 1; v <= m_nvoxels; v++)
    {
        if (m_continueFromFile != "")
        {
            m_fwd_post[v - 1] = resultMVNs.at(v - 1)->GetSubmatrix(1, m_num_params);

            assert(m_num_params + nNoiseParams == resultMVNs.at(v - 1)->GetSize());
            m_noise_post[v - 1] = noise->NewParams();
            m_noise_post[v - 1]->InputFromMVN(resultMVNs.at(v - 1)->GetSubmatrix(
                m_num_params + 1, m_num_params + nNoiseParams));
        }
        else
        {
            m_fwd_post[v - 1] = *initialFwdPosterior;
            m_noise_post[v - 1] = initialNoisePosterior->Clone();
        }

        if (m_locked_linear)
        {
            lockedLinearCentres.Column(v) = lockedLinearDists.at(v - 1)->means.Rows(1, m_num_params);
            m_lin_model[v - 1].ReCentre(lockedLinearCentres.Column(v));
        }
        else
        {
            m_lin_model[v - 1].ReCentre(m_fwd_post[v - 1].means);
        }

        if (m_save_without_prior)
        {
            m_fwd_post_no_prior.at(v - 1) = new MVNDist();
        }

        m_noise_prior[v - 1] = initialNoisePrior->Clone();
        noise->Precalculate(*m_noise_post[v - 1], *m_noise_prior[v - 1],
            m_origdata->Column(v));
    }
}

void SpatialVariationalBayes::SetupStSMatrix()
{
    assert((int)m_neighbours.size() == m_nvoxels);

    const double tiny = 1e-6;
    WARN_ONCE(
        "SpatialVariationalBayes::SetupStSMatrix Using 'S' prior with fast-calculation method and constant "
        "diagonal weight of "
        + stringify(tiny));
    // NEW METHOD
    {
        LOG << "SpatialVariationalBayes::SetupStSMatrix Attempting to allocate, m_nvoxels = " << m_nvoxels << endl;
        m_sts.ReSize(m_nvoxels);
        m_sts = 0;
        for (int v = 1; v <= m_nvoxels; v++)
        {
            // Number of neighbours v has
            int Nv = m_neighbours[v - 1].size();

            // Diagonal value = N + (N+tiny)^2
            m_sts(v, v) = Nv + (Nv + tiny) * (Nv + tiny);

            // Off-diagonal value = num 2nd-order neighbours (with duplicates) -
            // Aij(Ni+Nj+2*tiny)
            for (vector<int>::iterator nidIt = m_neighbours[v - 1].begin();
                 nidIt != m_neighbours[v - 1].end(); ++nidIt)
            {
                if (v < *nidIt)
                    m_sts(v, *nidIt) -= Nv + m_neighbours[*nidIt - 1].size() + 2 * tiny;
            }
            for (vector<int>::iterator nidIt = m_neighbours2[v - 1].begin();
                 nidIt != m_neighbours2[v - 1].end(); ++nidIt)
            {
                if (v < *nidIt)
                    m_sts(v, *nidIt) += 1;
            }
        }
        LOG << "SpatialVariationalBayes::Done generating StS matrix (New method)" << endl;
    }
}

void SpatialVariationalBayes::UpdateAkmean(DiagonalMatrix &akmean)
{
    LOG << "SpatialVariationalBayes::UpdateAkmean" << endl;
    assert(akmean.Nrows() == m_num_params);
    // Update spatial normalization term

    const double tiny = 0; // turns out to be no longer necessary.

    // Collect gk, wk, sigmak across all voxels
    DiagonalMatrix gk(m_num_params); // gk = 0.0/0.0;
    for (int k = 1; k <= m_num_params; k++)
    {
        ColumnVector wk(m_nvoxels);       // wk = 0.0/0.0;
        DiagonalMatrix sigmak(m_nvoxels); // sigmak=0.0/0.0;
        for (int v = 1; v <= m_nvoxels; v++)
        {
            wk(v) = m_fwd_post.at(v - 1).means(k);
            sigmak(v, v) = m_fwd_post.at(v - 1).GetCovariance()(k, k);
        }
        //    wk = wk(:);
        //    tmp1 = trace(diag(sigmak)*S'*S)

        // LOG << "m_fwd_post[0].GetCovariance():" <<
        // m_fwd_post[0].GetCovariance();

        // Update from Penny05:
        // 1/gk = 0.5*Trace[Sigmak*S'*S] + 0.5*wk'*S'*S*wk + 1/q1
        // hk = N/2 + q2
        // To do the MRF, just replace S'*S with S.

        switch (m_shrinkage_type)
        {
        case 'Z':
        {
            assert(m_save_without_prior);
            assert(m_sts.Nrows() == m_nvoxels);

            // Prior used by penny:
            //		  double q1 = 10, q2 = 1;
            // Noninformative prior:
            double q1 = 1e12, q2 = 1e-12;

            WARN_ONCE("SpatialVariationalBayes::UpdateAkmean Hyperpriors on S prior: using q1 == " + stringify(q1) + ", q2 == " + stringify(q2));

            gk(k) = 1 / (0.5 * (sigmak * m_sts).Trace() + (wk.t() * m_sts * wk).AsScalar() + 1 / q1);

            akmean(k) = gk(k) * (0.5 * m_nvoxels + q2);
        }
        // only output is akmean(k)
        break;

        case 'p':
        case 'P':
        case 'm':
        case 'M':
        case 'S':

        {
            // The following calculates Tr[Sigmak*S'*S]
            // using the fact that this == sum(diag(sigmak) .* diag(S'*S))
            // (since sigmak is diagonal!)

            double tmp1 = 0.0;
            for (int v = 1; v <= m_nvoxels; v++)
            {
                int nn = m_neighbours.at(v - 1).size();
                if (m_shrinkage_type == 'm') // useMRF)
                    tmp1 += sigmak(v, v) * m_spatial_dims * 2;
                else if (m_shrinkage_type == 'M') // useMRF2)
                    tmp1 += sigmak(v, v) * (nn + 1e-8);
                else if (m_shrinkage_type == 'p')
                    tmp1 += sigmak(v, v) * (4 * m_spatial_dims * m_spatial_dims + nn);
                else if (m_shrinkage_type == 'S')
                    tmp1 += sigmak(v, v) * ((nn + 1e-6) * (nn + 1e-6) + nn);
                else // P
                    tmp1 += sigmak(v, v) * ((nn + tiny) * (nn + tiny) + nn);
            }

            ColumnVector Swk = tiny * wk;

            if (m_shrinkage_type == 'S')
                Swk = 1e-6 * wk;

            for (int v = 1; v <= m_nvoxels; v++)
            {
                for (vector<int>::iterator v2It = m_neighbours[v - 1].begin();
                     v2It != m_neighbours.at(v - 1).end(); ++v2It)
                {
                    Swk(v) += wk(v) - wk(*v2It);
                }
                //		if (useDirichletBC || useMRF) // but not useMRF2
                if (m_shrinkage_type == 'p' || m_shrinkage_type == 'm')
                    Swk(v) += wk(v) * (m_spatial_dims * 2 - m_neighbours.at(v - 1).size());
                // Do nothing for 'S'
            }
            double tmp2 = Swk.SumSquare(); //(Swk.t() * Swk).AsScalar();

            //	    if (useMRF || useMRF2) // overwrite this for MRF
            if (m_shrinkage_type == 'm' || m_shrinkage_type == 'M')
                tmp2 = DotProduct(Swk, wk);

            LOG << "SpatialVariationalBayes::UpdateAkmean k=" << k << ", tmp1=" << tmp1 << ", tmp2=" << tmp2 << endl;
            // LOG << Swk.t();

            //  gk(k) = 1/(0.5*tmp1 + 0.5*tmp2 + 1/10)
            gk(k, k) = 1 / (0.5 * tmp1 + 0.5 * tmp2 + 0.1); // prior q1 == 10 (1/q1 == 0.1)
            //  end

            akmean(k) = gk(k) * (m_nvoxels * 0.5 + 1.0); // prior q2 == 1.0
        }

        break;

        default:
            assert(false);
            break;
        }
    }

    DiagonalMatrix akmeanMax = akmean * m_spatial_speed;

    for (int k = 1; k <= akmean.Nrows(); k++)
    {
        if (akmean(k, k) < 1e-50)
        {
            LOG << "SpatialVariationalBayes::UpdateAkmean akmean(" << k << ") was " << akmean(k, k) << endl;
            WARN_ONCE("SpatialVariationalBayes::UpdateAkmean akmean value was tiny!");
            akmean(k, k) = 1e-50; // prevent crashes
        }

        if (akmeanMax(k) < 0.5)
        {
            akmeanMax(k) = 0.5; // totally the wrong scaling.. oh well
        }

        if (m_spatial_speed > 0 && akmean(k) > akmeanMax(k))
        {
            LOG << "SpatialVariationalBayes::UpdateAkmean Rate-limiting the increase on akmean " << k << ": was "
                << akmean(k) << ", now " << akmeanMax(k) << endl;
            akmean(k) = akmeanMax(k);
        }
    }

    LOG << "SpatialVariationalBayes::UpdateAkmean New akmean: " << akmean.AsColumn().t();
}

void SpatialVariationalBayes::UpdateDeltaRho(DiagonalMatrix &delta, DiagonalMatrix &rho, const DiagonalMatrix &akmean, bool first_iter)
{
    for (int k = 1; k <= m_num_params; k++)
    {
        // if(k<6){LOG_ERR("Skipping parameter "<<k<<endl);continue;}
        LOG << "SpatialVariationalBayes::UpdateDeltaRho Optimizing for parameter " << k << endl;

        char type = m_prior_types[k - 1].m_type;

        // Each type should issue exactly one line to the logfile of the form
        // SpatialPrior on k type ? ?? : x y z
        // ? = single-character type
        // ?? = any other subtype info (optional, free-form but no : character)
        // x, y, z = numerical parameters (e.g. delta, rho, 0)
        switch (type)
        {
        case '-':
        case 'N':
        case 'I':
        case 'A':
            // Nonspatial priors
            delta(k) = 0;
            rho(k) = 0; // no shrinkage/scaling factor either
            LOG << "SpatialVariationalBayes::SpatialPrior " << k << " type " << type << endl;
            break;

        //	  case 'F':
        //	    delta(k) = m_fixed_delta;
        //	    rho(k) = m_fixed_rho;
        //	    break;

        case 'm':
        case 'M':
        case 'p':
        case 'P':
        case 'S':

            assert(type == m_shrinkage_type);
            // fill with invalid values:
            delta(k) = -3;
            rho(k) = 1234.5678;
            LOG << "SpatialVariationalBayes::SpatialPrior " << k << " type " << type << endl;

            break;

        default:
            throw InvalidOptionValue("param-spatial-priors", stringify(type),
                "Invalid spatial prior type");

        case 'R':
        case 'D':
        case 'F':
            // Reorganize data by parameter (rather than by voxel)
            DiagonalMatrix covRatio(m_nvoxels);
            ColumnVector meanDiffRatio(m_nvoxels);
            const double priorCov = initialFwdPrior->GetCovariance()(k, k);
            const double priorCovSqrt = sqrt(priorCov);
            const double priorMean = initialFwdPrior->means(k);

            for (int v = 1; v <= m_nvoxels; v++)
            {
                // Isolate just the dimensionless quantities we need

                // Penny:
                covRatio(v, v) = m_fwd_post.at(v - 1).GetCovariance()(k, k) / priorCov;
                // Hacky:
                // LOG_ERR("WARNING: Using hacky covRatio calculation (precision
                // rather than covariance!\n!");
                //		covRatio(v,v) = 1 /
                // m_fwd_post.at(v-1).GetPrecisions()(k,k) / priorCov;

                meanDiffRatio(v) = (m_fwd_post.at(v - 1).means(k) - priorMean) / priorCovSqrt;
            }

            //	SymmetricMatrix covRatioSupplemented(m_nvoxels);
            // Recover the off-diagonal elements of m_fwd_prior/priorCov

            /* TURN COVRATIOSUPPLEMENTED BACK INTO COVRATIO!

                    if (Sinvs.at(k-1).Nrows() != 0)
                    {
                    // Implicitly use all zeroes for first iteration
                    assert(Sinvs.at(k-1).Nrows() == m_nvoxels);
                    covRatioSupplemented = Sinvs.at(k-1);

                    for (int v=1; v<=m_nvoxels; v++)
                    {
                    //		LOG << "Should be one: " <<
    initialFwdPrior->GetPrecisions()(k,k) * Sinvs.at(k-1)(v,v) /
    m_fwd_prior.at(v-1).GetPrecisions()(k,k) << endl;
                    assert(initialFwdPrior->GetPrecisions()(k,k)*Sinvs.at(k-1)(v,v)
    == m_fwd_prior.at(v-1).GetPrecisions()(k,k));

                    covRatioSupplemented(v,v) = 0;
                    }
                    }

                    covRatioSupplemented += covRatio.i();
                    covRatioSupplemented = covRatioSupplemented.i();
                    */

            // INSTEAD DO THIS:
            //	covRatioSupplemented = covRatio;
            if (first_iter && !m_update_first_iter)
            {
                if (type == 'F' && m_brute_force_delta_search)
                    LOG << "SpatialVariationalBayes::UpdateDeltaRho Doing calc on first iteration, just because it's F and "
                           "bruteForceDeltaSearch is on.  Temporary hack!\n";
                else
                    break; // skip the updates
            }

            DiagonalMatrix deltaMax = delta * m_spatial_speed;

            if (type == 'R')
            {
                if (m_always_inital_delta_guess > 0)
                    delta(k) = m_always_inital_delta_guess;
                if (m_use_evidence)
                {
                    WARN_ALWAYS("SpatialVariationalBayes::UpdateDeltaRho Using R... mistake??");
                    delta(k) = OptimizeEvidence(m_fwd_post_no_prior, k,
                        initialFwdPrior.get(), delta(k), true,
                        &rho(k));
                    LOG << "SpatialVariationalBayes::SpatialPrior " << k << " type R eo : " << delta(k)
                        << " " << rho(k) << " 0" << endl;
                }
                else
                {
                    WARN_ALWAYS("SpatialVariationalBayes::UpdateDeltaRho Using R without EO... mistake??");
                    // Spatial priors with rho & delta
                    delta(k) = OptimizeSmoothingScale(covRatio, meanDiffRatio,
                        delta(k), &rho(k), true);
                    LOG << "SpatialVariationalBayes::SpatialPrior " << k << " type R vb : " << delta(k)
                        << " " << rho(k) << " 0" << endl;
                }
            }
            else if (type == 'D')
            {
                if (m_always_inital_delta_guess > 0)
                    delta(k) = m_always_inital_delta_guess;

                // Spatial priors with only delta
                if (m_use_evidence)
                {
                    delta(k) = OptimizeEvidence(m_fwd_post_no_prior, k,
                        initialFwdPrior.get(), delta(k));
                    LOG << "SpatialVariationalBayes::SpatialPrior " << k << " type D eo : " << delta(k)
                        << " 0 0" << endl;
                }
                else
                {
                    WARN_ALWAYS("SpatialVariationalBayes::UpdateDeltaRho Using D without EO... mistake??");
                    delta(k) = OptimizeSmoothingScale(covRatio, meanDiffRatio,
                        delta(k), &rho(k), false);
                    LOG << "SpatialVariationalBayes::SpatialPrior " << k << " type D vb : " << delta(k)
                        << " 0 0" << endl;
                }
            }
            else // type == 'F'
            {
                delta(k) = m_fixed_delta;
                rho(k) = m_fixed_rho;

                // The following does nothing BUT it's neccessary to
                // make the bruteForceDeltaEstimates work.
                double newDelta = OptimizeSmoothingScale(
                    covRatio, meanDiffRatio, delta(k), &rho(k), false, false);
                assert(newDelta == m_fixed_delta);
                assert(rho(k) == m_fixed_rho);
                deltaMax(k) = delta(k);
                LOG << "SpatialVariationalBayes::SpatialPrior " << k << " type F : " << delta(k) << " "
                    << rho(k) << " 0" << endl;
            }

            // enforce maxPrecisionIncreasePerIteration
            if (deltaMax(k) < 0.5)
                deltaMax(k) = 0.5;
            if (m_spatial_speed > 0 && delta(k) > deltaMax(k))
            {
                LOG << "SpatialVariationalBayes::Rate-limiting the increase on delta " << k << ": was "
                    << delta(k);
                delta(k) = deltaMax(k);
                LOG << ", now " << delta(k) << endl;

                // Re-evaluate rho, for this delta
                double newDelta = OptimizeSmoothingScale(
                    covRatio, meanDiffRatio, delta(k), &rho(k), type == 'R', false);
                assert(newDelta == delta(k)); // a quick check
            }
            // default: dealt with earlier.
        }

        LOG << "SpatialVariationalBayes::UpdateDeltaRho delta(k) = " << delta(k) << ", rho(k) == " << rho(k)
            << endl;
    }
}

// CALCULATE THE C^-1 FOR THE NEW DELTAS
void SpatialVariationalBayes::CalculateCinv(vector<SymmetricMatrix> &Sinvs, DiagonalMatrix &delta, DiagonalMatrix &rho, DiagonalMatrix &akmean)
{
    LOG << "SpatialVariationalBayes::CalculateCinv" << endl;
    for (int k = 1; k <= m_num_params; k++)
    {
        if (delta(k) >= 0)
        {
            Sinvs.at(k - 1) = m_covar.GetCinv(delta(k)) * exp(rho(k));

            if (delta(k) == 0 && m_save_without_prior)
            {
                assert((m_prior_types[k - 1].m_type == 'N') | (m_prior_types[k - 1].m_type == 'I') | (m_prior_types[k - 1].m_type == 'A'));
                Sinvs.at(k - 1) = IdentityMatrix(m_nvoxels);
            }

            assert(SP(initialFwdPrior->GetPrecisions(),
                       IdentityMatrix(m_num_params) - 1)
                       .MaximumAbsoluteValue()
                == 0);
            Sinvs[k - 1] *= initialFwdPrior->GetPrecisions()(k, k);
        }

        if (delta(k) < 0 && m_save_without_prior)
        {
            assert(m_prior_types[k - 1].m_type == m_shrinkage_type);

            if (m_shrinkage_type == 'S')
            {
                assert(m_sts.Nrows() == m_nvoxels);
                Sinvs.at(k - 1) = m_sts * akmean(k);
            }
            else
            {
                assert(m_shrinkage_type == 'p');

                // Build up the second-order matrix directly, row-by-row
                Matrix sTmp(m_nvoxels, m_nvoxels);
                sTmp = 0;

                for (int v = 1; v <= m_nvoxels; v++)
                {
                    // self = (2*Ndim)^2 + (nn)
                    sTmp(v, v) = 4 * m_spatial_dims * m_spatial_dims; // nn added later

                    // neighbours = (2*Ndim) * -2
                    for (vector<int>::iterator nidIt = m_neighbours[v - 1].begin();
                         nidIt != m_neighbours[v - 1].end(); ++nidIt)
                    {
                        int nid = *nidIt; // neighbour ID (voxel number)
                        assert(sTmp(v, nid) == 0);
                        sTmp(v, nid) = -2 * 2 * m_spatial_dims;
                        sTmp(v, v) += 1;
                    }

                    // neighbours2 = 1 (for each appearance)
                    for (vector<int>::iterator nidIt = m_neighbours2[v - 1].begin();
                         nidIt != m_neighbours2[v - 1].end(); ++nidIt)
                    {
                        int nid2 = *nidIt;  // neighbour ID (voxel number)
                        sTmp(v, nid2) += 1; // not =1, because duplicates are ok.
                    }
                }

                // Store back into Sinvs[k-1] and apply akmean(k)
                assert(sTmp == sTmp.t());
                Sinvs.at(k - 1) << sTmp;
                Sinvs[k - 1] *= akmean(k);
            }
        }
    }
}

void SpatialVariationalBayes::DoSimEvidence(vector<SymmetricMatrix> &Sinvs)
{
    LOG << "SpatialVariationalBayes::DoSimEvidence";

    // Re-estimate m_fwd_prior for all voxels simultaneously,
    // based on the full covariance matrix

    // Check it's the simple case (haven't coded up the correction
    // factors yet)
    //	assert(initialFwdPrior->GetPrecisions() ==
    // IdentityMatrix(m_num_params)); // now part of Sinvs
    //	assert(initialFwdPrior->means == -initialFwdPrior->means);  //
    // but this still applies

    if (!(initialFwdPrior->means == -initialFwdPrior->means))
        WARN_ALWAYS("SpatialVariationalBayes::DoSimEvidence Quick hack to avoid assertion with initialFwdPrior->means != 0");

    SymmetricMatrix SigmaInv(m_num_params * m_nvoxels);
    //	SymmetricMatrix Sigma(m_num_params*m_nvoxels);
    ColumnVector Mu(m_num_params * m_nvoxels);

    // These matrices consist of NxN matrices blocked together
    // so parameter k, voxel v is in row (or col): v + (k-1)*m_num_params
    SymmetricMatrix Ci = -999 * IdentityMatrix(m_nvoxels * m_num_params);
    SymmetricMatrix XXtr = -999 * IdentityMatrix(m_nvoxels * m_num_params);
    ColumnVector XYtr(m_nvoxels * m_num_params);
    XYtr = -999;

    // Build Ci
    for (int k = 1; k <= m_num_params; k++)
    {
        Ci.SymSubMatrix(m_nvoxels * (k - 1) + 1, m_nvoxels * k) = Sinvs[k - 1];
        // off-diagonal blocks are zero, by definition of the our priors
        // (priors between parameters are independent)
    }

    // Build XXtr and XYtr
    for (int v = 1; v <= m_nvoxels; v++)
    {
        const SymmetricMatrix &tmp = m_fwd_post_no_prior[v - 1]->GetPrecisions();
        ColumnVector tmp2 = tmp * (m_fwd_post_no_prior[v - 1]->means - initialFwdPrior->means);
        for (int k1 = 1; k1 <= m_num_params; k1++)
        {
            XYtr(v + (k1 - 1) * m_nvoxels) = tmp2(k1);
            for (int k2 = 1; k2 <= m_num_params; k2++)
            {
                XXtr(v + (k1 - 1) * m_nvoxels, v + (k2 - 1) * m_nvoxels) = tmp(k1, k2);
            }
        }
        // XXtr != 0 only where the row and col refer to the same voxel.
    }

    SigmaInv = XXtr + Ci;
    Mu = SigmaInv.i() * XYtr;

    for (int v = 1; v <= m_nvoxels; v++)
    {
        ColumnVector muBefore = m_fwd_post[v - 1].means - initialFwdPrior->means;

        assert(m_full_eo_first_param == 1);

        for (int k = 1; k <= m_num_params; k++)
            m_fwd_post[v - 1].means(k) = Mu(v + (k - 1) * m_nvoxels) + initialFwdPrior->means(k);

        //	    if ((muBefore -
        // m_fwd_post[v-1].means).MaximumAbsoluteValue() > 1e-10)
        //	      LOG << "mBef = " << muBefore.t() << "mAft = " <<
        // m_fwd_post[v-1].means.t();

        if (m_use_covar_marginals_not_precisions)
        {
            SymmetricMatrix Sigma = SigmaInv.i();

            SymmetricMatrix cov = m_fwd_post[v - 1].GetCovariance();
            SymmetricMatrix covOld = cov;

            WARN_ONCE("SpatialVariationalBayes::DoSimEvidence Full simultaneous diagonal thingy -- now in covariances!");

            for (int k1 = 1; k1 <= m_num_params; k1++)
            {
                for (int k2 = 1; k2 <= m_num_params; k2++)
                {
                    cov(k1, k2) = Sigma(v + (k1 - 1) * m_nvoxels, v + (k2 - 1) * m_nvoxels);
                }
            }
            if ((cov - covOld).MaximumAbsoluteValue() > 1e-10)
                LOG << "SpatialVariationalBayes::DoSimEvidence covBefore: " << covOld.AsColumn().t()
                    << "covAfter: " << cov.AsColumn().t();
            m_fwd_post[v - 1].SetCovariance(cov);
        }
        else
        {
            SymmetricMatrix prec = m_fwd_post[v - 1].GetPrecisions();

            SymmetricMatrix precOld = prec;
            // LOG << "precBefore:\n" << prec;

            WARN_ONCE("SpatialVariationalBayes::DoSimEvidence Full simultaneous diagonal thingy");
            // LOG << "prec = \n" << prec << endl;
            for (int k1 = 1; k1 <= m_num_params; k1++)
            {
                for (int k2 = 1; k2 <= m_num_params; k2++)
                {
                    prec(k1, k2) = SigmaInv(v + (k1 - 1) * m_nvoxels, v + (k2 - 1) * m_nvoxels);
                }
            }

            if ((prec - precOld).MaximumAbsoluteValue() > 1e-10)
                LOG << "SpatialVariationalBayes::DoSimEvidence precBefore: " << precOld.AsColumn().t()
                    << "precAfter: " << prec.AsColumn().t();

            // if ((m_fwd_post[v-1].GetPrecisions() -
            // prec).MaximumAbsoluteValue() > 1e-10)
            //  LOG << "pbef: " << m_fwd_post[v-1].GetPrecisions() << "paft:
            //  " << prec;

            m_fwd_post[v - 1].SetPrecisions(prec);
        }
    }
}

void SpatialVariationalBayes::DoFullEvidence(vector<SymmetricMatrix> &Sinvs)
{
    //	assert(!m_use_covar_marginals_not_precisions);
    // Covariance marginals are broken below, and I think they're
    // rubbish anyway
    LOG << "SpatialVariationalBayes::DoFullEvidence using " + string(m_use_covar_marginals_not_precisions
                                                                      ? "covariances."
                                                                      : "precisions.");

    // Re-estimate m_fwd_prior for all voxels simultaneously,
    // based on the full covariance matrix

    // Check it's the simple case (haven't coded up the correction
    // factors yet)
    //	assert(initialFwdPrior->GetPrecisions() ==
    // IdentityMatrix(m_num_params)); // now part of Sinvs
    //	assert(initialFwdPrior->means == -initialFwdPrior->means);  //
    // but this still applies

    vector<SymmetricMatrix> SigmaInv(m_num_params);
    vector<SymmetricMatrix> Sigma(m_num_params);
    vector<ColumnVector> Mu(m_num_params);

    for (int k = 1; k <= m_num_params; k++)
    {
        const SymmetricMatrix &Ci = Sinvs[k - 1];
        SymmetricMatrix XXtr(m_nvoxels);
        ColumnVector XYtr(m_nvoxels);

        ColumnVector XXtrMuOthers(m_nvoxels);

        // Initialize to junk values
        XXtr = -999 * IdentityMatrix(m_nvoxels);
        XYtr = -999;

        for (int v = 1; v <= m_nvoxels; v++)
        {
            const SymmetricMatrix &tmp = m_fwd_post_no_prior[v - 1]->GetPrecisions();
            XXtr(v, v) = tmp(k, k);

            ColumnVector tmp2 = tmp * (m_fwd_post_no_prior[v - 1]->means - initialFwdPrior->means);
            XYtr(v) = tmp2(k);

            //		ColumnVector MuOthers = m_fwd_post[v-1].means;
            ColumnVector MuOthers = m_fwd_post[v - 1].means - initialFwdPrior->means;
            MuOthers(k) = 0;
            ColumnVector tmp3 = tmp * MuOthers;
            XXtrMuOthers(v) = tmp3(k);

            WARN_ONCE(
                "SpatialVariationalBayes::DoFullEvidence Corrected mistake in useFullEvidenceOptimization: "
                "initialFwdPrior->means (not k)");
            // Also notice the subtle difference above: MuOthers uses the actual
            // posterior means, while XYtr uses
            // the priorless posterior means.
            // Also, the above XXtr do NOT include the correction for non-N(0,1)
            // initialFwdPriors, because the
            // Sinvs (and hence Ci etc) already include this correction.  This is
            // different from the DerivEdDelta
            // calculation which uses Cs with 1 on the diagonal (originally to
            // facilitate reuse across different rho
            // values, now just confusing).
        }

        //	    ColumnVector tmp4(m_nvoxels);
        //	    tmp4 = initialFwdPrior->means(k);
        //	    ColumnVector CiMu0 = Ci * tmp4;

        {
            SigmaInv.at(k - 1) = XXtr + Ci;
        }
        {
            Sigma.at(k - 1) = SigmaInv[k - 1].i();
        }
        {
            //	      Mu.at(k-1) = Sigma[k-1] * (XYtr - XXtrMuOthers);
            //	      Mu.at(k-1) = Sigma[k-1] * (XYtr - XXtrMuOthers + CiMu0);
            Mu.at(k - 1) = Sigma[k - 1] * (XYtr - XXtrMuOthers);
        }
    }
    for (int v = 1; v <= m_nvoxels; v++)
    {
        ColumnVector muBefore = m_fwd_post[v - 1].means;

        for (int k = m_full_eo_first_param; k <= m_num_params; k++)
            m_fwd_post[v - 1].means(k) = Mu[k - 1](v) + initialFwdPrior->means(k);

        //	    if ((muBefore -
        // m_fwd_post[v-1].means).MaximumAbsoluteValue() > 1e-10)
        //	      {
        //		LOG << "mBef = " << muBefore.t() << "mAft = " <<
        // m_fwd_post[v-1].means.t();
        //	      }

        if (m_use_covar_marginals_not_precisions)
        {
            SymmetricMatrix cov = SP(m_fwd_post[v - 1].GetCovariance(),
                IdentityMatrix(m_num_params));
            WARN_ONCE("SpatialVariationalBayes::DoFullEvidence Covariance diagonal thingy");
            // LOG << "cov = \n" << cov << endl;

            for (int k = m_full_eo_first_param; k <= m_num_params; k++)
                cov(k, k) = Sigma[k - 1](v, v);

            m_fwd_post[v - 1].SetCovariance(cov);
        }
        else if (m_keep_param_covars)
        {
            WARN_ONCE("SpatialVariationalBayes::DoFullEvidence Keeping inter-parameter covariances from VB!");
        }
        else
        {
            SymmetricMatrix prec = SP(m_fwd_post[v - 1].GetPrecisions(),
                IdentityMatrix(m_num_params));

            SymmetricMatrix precOld = prec;
            // LOG << "precBefore:\n" << prec;

            WARN_ONCE("SpatialVariationalBayes::DoFullEvidence Precision diagonal thingy");
            // LOG << "prec = \n" << prec << endl;
            for (int k = m_full_eo_first_param; k <= m_num_params; k++)
                prec(k, k) = SigmaInv[k - 1](v, v);

            if ((prec - precOld).MaximumAbsoluteValue() > 1e-10)
                LOG << "SpatialVariationalBayes::DoFullEvidence precBefore: " << precOld.AsColumn().t()
                    << "precAfter: " << prec.AsColumn().t();

            // if ((m_fwd_post[v-1].GetPrecisions() -
            // prec).MaximumAbsoluteValue() > 1e-10)
            //  LOG << "pbef: " << m_fwd_post[v-1].GetPrecisions() << "paft:
            //  " << prec;

            m_fwd_post[v - 1].SetPrecisions(prec);
        }
        assert(m_fwd_post[v - 1].GetSize() == m_num_params);
    }
}

void SpatialVariationalBayes::SetFwdPriorShrinkageTypeS(int v, const NEWMAT::DiagonalMatrix &akmean)
{
    WARN_ONCE("SpatialVariationalBayes::SetFwdPriorShrinkageTypeS Using new S VB spatial thingy");

    assert(m_sts.Nrows() == m_nvoxels);

    double weight = 1e-6; // weakly pulled to zero
    ColumnVector contrib(m_num_params);
    contrib = 0;

    for (int i = 1; i <= m_nvoxels; i++)
    {
        if (v != i)
        {
            weight += m_sts(v, i);
            contrib += m_sts(v, i) * m_fwd_post[i - 1].means;
        }
    }

    DiagonalMatrix spatialPrecisions;
    spatialPrecisions = akmean * m_sts(v, v);

    m_fwd_prior[v - 1].SetPrecisions(spatialPrecisions);
    m_fwd_prior[v - 1].means = contrib / weight;
}

void SpatialVariationalBayes::SetFwdPriorShrinkageType(int v, const NEWMAT::DiagonalMatrix &akmean)
{
    const double tiny = 0; // turns out to be no longer necessary.
    double weight8 = 0;    // weighted +8
    ColumnVector contrib8(m_num_params);
    contrib8 = 0.0;
    for (vector<int>::iterator nidIt = m_neighbours[v - 1].begin();
         nidIt != m_neighbours[v - 1].end(); ++nidIt)
    // iterate over neighbour ids
    {
        int nid = *nidIt;
        const MVNDist &neighbourPost = m_fwd_post[nid - 1];
        contrib8 += 8 * neighbourPost.means;
        weight8 += 8;
    }

    double weight12 = 0; // weighted -1, may be duplicated
    ColumnVector contrib12(m_num_params);
    contrib12 = 0.0;
    for (vector<int>::iterator nidIt = m_neighbours2[v - 1].begin();
         nidIt != m_neighbours2[v - 1].end(); ++nidIt)
    // iterate over neighbour ids
    {
        int nid = *nidIt;
        const MVNDist &neighbourPost = m_fwd_post[nid - 1];
        contrib12 += -neighbourPost.means;
        weight12 += -1;
    }

    // Set prior mean & precisions
    int nn = m_neighbours[v - 1].size();

    if (m_shrinkage_type == 'p')
    {
        //	LOG << nn << " -> " << 2*spatialDims << endl;
        assert(nn <= m_spatial_dims * 2);
        // nn = spatialDims*2;
        weight8 = 8 * 2 * m_spatial_dims;
        weight12 = -1 * (4 * m_spatial_dims * m_spatial_dims - nn);
    }

    DiagonalMatrix spatialPrecisions;

    if (m_shrinkage_type == 'P')
        spatialPrecisions = akmean * ((nn + tiny) * (nn + tiny) + nn);
    else if (m_shrinkage_type == 'm')
        spatialPrecisions = akmean * m_spatial_dims * 2;
    else if (m_shrinkage_type == 'M')
        spatialPrecisions = akmean * (nn + 1e-8);
    else if (m_shrinkage_type == 'p')
        spatialPrecisions = akmean * (4 * m_spatial_dims * m_spatial_dims + nn);
    else if (m_shrinkage_type == 'S')
    {
        spatialPrecisions = akmean * ((nn + 1e-6) * (nn + 1e-6) + nn);
        WARN_ONCE("SpatialVariationalBayes::SetFwdPriorShrinkageType Using a hacked-together VB version of the 'S' prior");
    }

    //	    if (useDirichletBC || useMRF)
    if (m_shrinkage_type == 'p' || m_shrinkage_type == 'm')
    {
        //	LOG_ERR("Penny-style DirichletBC priors -- ignoring
        // initialFwdPrior completely!\n");
        m_fwd_prior[v - 1].SetPrecisions(spatialPrecisions);
    }
    else
    {
        m_fwd_prior[v - 1].SetPrecisions(initialFwdPrior->GetPrecisions() + spatialPrecisions);
    }

    ColumnVector mTmp(m_num_params);

    if (weight8 != 0)
        mTmp = (contrib8 + contrib12) / (weight8 + weight12);
    else
        mTmp = 0;

    if (m_shrinkage_type == 'm')                    // useMRF) // overwrite this for MRF
        mTmp = contrib8 / (8 * m_spatial_dims * 2); // note: Dirichlet BCs on MRF
    if (m_shrinkage_type == 'M')                    // useMRF2)
        mTmp = contrib8 / (8 * (nn + 1e-8));

    // equivalent, when non-spatial priors are very weak:
    //    m_fwd_prior[v-1].means = mTmp;

    m_fwd_prior[v - 1].means = m_fwd_prior[v - 1].GetCovariance() * (spatialPrecisions * mTmp + initialFwdPrior->GetPrecisions() * initialFwdPrior->means);

    if (m_shrinkage_type == 'm' || m_shrinkage_type == 'M')
        m_fwd_prior[v - 1].means = m_fwd_prior[v - 1].GetCovariance() * spatialPrecisions * mTmp; // = mTmp;
}

double SpatialVariationalBayes::SetFwdPrior(int v, const vector<SymmetricMatrix> &Sinvs, bool isFirstIteration)
{
    // Use the new spatial priors
    // Marginalize out all the other voxels
    double Fard = 0;

    DiagonalMatrix spatialPrecisions(m_num_params);
    ColumnVector weightedMeans(m_num_params);

    ColumnVector priorMeans(m_num_params);
    // default is to get these from intialFwdPrior
    // this is overwritten for I priors
    // or ignored for spatial priors
    priorMeans = initialFwdPrior->means;

    for (int k = 1; k <= m_num_params; k++)
    {
        if (m_prior_types[k - 1].m_type == m_shrinkage_type)
        {
            spatialPrecisions(k) = -9999;
            weightedMeans(k) = -9999;
            continue;
        }
        else if (m_prior_types[k - 1].m_type == 'A')
        {
            if (isFirstIteration)
            {
                spatialPrecisions(k) = initialFwdPrior->GetPrecisions()(k, k);
                weightedMeans(k) = initialFwdPrior->means(k);
                // Fard = 0;
            }
            else
            {
                double ARDparam = 1 / m_fwd_post[v - 1].GetPrecisions()(k, k) + m_fwd_post[v - 1].means(k) * m_fwd_post[v - 1].means(k);
                spatialPrecisions(k) = 1 / ARDparam;
                weightedMeans(k) = 0;
                Fard -= 2.0 * log(2.0 / ARDparam);
            }
            continue;
        }
        else if (m_prior_types[k - 1].m_type == 'N')
        {
            // special case because Sinvs is 0x0, but should actually
            // be the identity matrix.
            spatialPrecisions(k) = initialFwdPrior->GetPrecisions()(k, k);
            assert(SP(initialFwdPrior->GetPrecisions(),
                       IdentityMatrix(m_num_params) - 1)
                       .MaximumAbsoluteValue()
                == 0);

            // Don't worry, this is multiplied by initialFwdPrior later
            weightedMeans(k) = 0;
            continue;
        }
        else if (m_prior_types[k - 1].m_type == 'I')
        {
            // get means from image prior MVN
            priorMeans(k) = m_prior_types[k - 1].m_image(v);

            // precisions in same way as 'N' prior (for time being!)
            spatialPrecisions(k) = initialFwdPrior->GetPrecisions()(k, k);
            assert(SP(initialFwdPrior->GetPrecisions(),
                       IdentityMatrix(m_num_params) - 1)
                       .MaximumAbsoluteValue()
                == 0);

            weightedMeans(k) = 0;
            continue;
        }

        spatialPrecisions(k) = Sinvs[k - 1](v, v);
        //	  double testWeights = 0;
        weightedMeans(k) = 0;
        for (int n = 1; n <= m_nvoxels; n++)
            if (n != v)
            {
                weightedMeans(k) += Sinvs[k - 1](n, v) * (m_fwd_post[n - 1].means(k) - initialFwdPrior->means(k));
                //	      testWeights += Cinvs[k-1](n,v);
            }
        //	  LOG_ERR("Parameter " << k << ", testWeights == " <<
        // testWeights << ", spatialPrecisions(k) == " << spatialPrecisions(k)
        //<< ", delta(k) == " << delta(k) << ", test2 == " << test2 << endl);
    }
    //      LOG_ERR("--------- end of voxel " << v << endl);

    assert(initialFwdPrior->GetPrecisions().Nrows() == spatialPrecisions.Nrows());
    // Should check that earlier!  It's possible for basis=1 and priors=2x2
    // to slip through.  TODO

    // Should check that initialFwdPrior was already diagonal --
    // this will cause real problems if it isn't!!
    // (Safe way: SP of the covariance matrices -- that'd force
    // diagonality while preserving individual variance.)

    // LOG << "Spatial precisions: " << spatialPrecisions;

    DiagonalMatrix finalPrecisions = spatialPrecisions;
    //	    SP(initialFwdPrior->GetPrecisions(),spatialPrecisions);

    // LOG << "initialFwdPrior->GetPrecisions() == " <<
    // initialFwdPrior->GetPrecisions();
    // LOG << "initialFwdPrior->GetCovariance() == " <<
    // initialFwdPrior->GetCovariance();

    ColumnVector finalMeans = priorMeans - spatialPrecisions.i() * weightedMeans;

    // LOG << "Final means and precisions: " << finalMeans <<
    // finalPrecisions;

    // Preserve the m_shrinkage_type ones from before.
    // They'd better be diagonal!
    for (int k = 1; k <= m_num_params; k++)
        if (m_prior_types[k - 1].m_type == m_shrinkage_type)
        {
            finalPrecisions(k) = m_fwd_prior[v - 1].GetPrecisions()(k, k);
            finalMeans(k) = m_fwd_prior[v - 1].means(k);
        }

    m_fwd_prior[v - 1].SetPrecisions(finalPrecisions);
    m_fwd_prior[v - 1].means = finalMeans;
    // Definitely a minus here.

    return Fard;
}

void SpatialVariationalBayes::DoCalculations(FabberRunData &allData)
{
    // extract data (and the coords) from allData for the (first) VB run
    // Rows are volumes
    // Columns are (time) series
    // num Rows is size of (time) series
    // num Cols is size of volumes
    m_origdata = &allData.GetMainVoxelData();
    m_coords = &allData.GetVoxelCoords();
    m_suppdata = &allData.GetVoxelSuppData();
    m_nvoxels = m_origdata->Ncols();
    int maxits = convertTo<int>(allData.GetStringDefault("max-iterations", "10"));

    // pass in some (dummy) data/coords here just in case the model relies upon it
    // use the first voxel values as our dummies FIXME this shouldn't really be
    // necessary, need to find way for model to know about the data beforehand.
    if (m_nvoxels > 0)
        PassModelData(1);

    // Only call DoCalculations once
    assert(resultMVNs.empty());
    assert(resultFs.empty());
    assert(resultMVNsWithoutPrior.empty());

    // Initialization:

    // Make the neighbours[] lists if required
    if (m_prior_types_str.find_first_of("mMpPSZ") != string::npos)
    {
        CalcNeighbours(allData.GetVoxelCoords());
    }

    // Make distance matrix if required
    if (m_prior_types_str.find_first_of("RDF") != string::npos)
    {
        // Note: really ought to know the voxel dimensions and multiply by those,
        // because CalcDistances expects an input in mm, not index.
        m_covar.CalcDistances(allData.GetVoxelCoords(), m_dist_measure);
    }

    SetupPerVoxelDists(allData);
    if (allData.GetBool("output-only"))
    {
        // Do no calculations - now we have set resultMVNs we can finish
        return;
    }

    // Make the spatial normalization parameters
    // akmean = 0*distsMaster.theta.means + 1e-8;
    DiagonalMatrix akmean(m_num_params);
    DiagonalMatrix delta(m_num_params);
    DiagonalMatrix rho(m_num_params);

    akmean = 1e-8;
    delta = m_fixed_delta; // Hard-coded initial value (in mm!)
    rho = 0;
    LOG << "SpatialVariationalBayes::Using initial value for all deltas: " << delta(1) << endl;

    vector<SymmetricMatrix> Sinvs(m_num_params);

    // FIXME can't calculate free energy with spatial VB yet 
    // This value never changes
    const double globalF = 1234.5678;

    // Cache for m_sts matrix in 'S' or 'Z' mode
    if ((m_shrinkage_type == 'S') || (m_shrinkage_type == 'Z'))
        SetupStSMatrix();

    m_conv->Reset();
    
    // MAIN ITERATION LOOP
    int it = 0;
    do
    {
        // Give an indication of the progress through the voxels;
        allData.Progress(it, maxits);

        // UPDATE SPATIAL SHRINKAGE PRIOR PARAMETERS
        if (m_shrinkage_type != '-' && ((it > 0) || m_update_first_iter))
        {
            UpdateAkmean(akmean);
        }
        
        UpdateDeltaRho(delta, rho, akmean, (it == 0));
        CalculateCinv(Sinvs, delta, rho, akmean);

        // ITERATE OVER VOXELS
        for (int v = 1; v <= m_nvoxels; v++)
        {
            PassModelData(v);

            if (m_continueFromFile == "")
            {
                // voxelwise initialisation - only if we dont have initial values from a
                // pre loaded MVN
                m_model->InitParams(m_fwd_post[v - 1]);
            }

            // Spatial VB mucks around with the fwd priors. Presumably to
            // incorporate the prior assumption of spatial uniformity

            // from simple_do_vb_ar1c_spatial.m
            // Note: this sets the priors as if all parameters were m_shrinkage_type.
            // We overwrite the non-m_shrinkage_type parameter priors later.
            if (m_shrinkage_type == 'S')
            {
                SetFwdPriorShrinkageTypeS(v, akmean);
            }
            else if (m_shrinkage_type != '-')
            {
                SetFwdPriorShrinkageType(v, akmean);
            }
            double Fard = SetFwdPrior(v, Sinvs, (it == 0));

            // The steps below are essentially the same as regular VB, although
            // the code looks different as the per-voxel dists are set up at the
            // start rather than as we go

            // Reference to free energy output so can be modified below
            double &F = resultFs.at(v - 1);
            if (m_needF)
            {
                F = noise->CalcFreeEnergy(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                    m_fwd_post[v - 1], m_fwd_prior[v - 1],
                    m_lin_model[v - 1], m_origdata->Column(v));
                F += Fard;
            }
            if (m_printF)
            {
                LOG << "SpatialVbInferenceTechnique::Fbefore == " << F << endl;
            }

            noise->UpdateTheta(*m_noise_post[v - 1], m_fwd_post[v - 1],
                m_fwd_prior[v - 1], m_lin_model[v - 1],
                m_origdata->Column(v),
                m_fwd_post_no_prior.at(v - 1));

            if (m_needF)
            {
                F = noise->CalcFreeEnergy(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                    m_fwd_post[v - 1], m_fwd_prior[v - 1],
                    m_lin_model[v - 1], m_origdata->Column(v));
                F += Fard;
                // Fard does NOT change because we haven't updated m_fwd_prior yet.
            }

            if (m_printF)
                LOG << "SpatialVbInferenceTechnique::Ftheta == " << F << endl;
        }

        // QUICK INTERRUPTION: Voxelwise calculations continue below.
        if (m_use_sim_evidence)
        {
            DoSimEvidence(Sinvs);
        }
        else if (m_use_full_evidence)
        {
            DoFullEvidence(Sinvs);
        }

        // Back to your regularly-scheduled voxelwise calculations
        for (int v = 1; v <= m_nvoxels; v++)
        {
            PassModelData(v);

            // Reference to free energy output so can be modified below
            double &F = resultFs.at(v - 1);

            noise->UpdateNoise(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                m_fwd_post[v - 1], m_lin_model[v - 1],
                m_origdata->Column(v));

            if (m_needF)
                F = noise->CalcFreeEnergy(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                    m_fwd_post[v - 1], m_fwd_prior[v - 1],
                    m_lin_model[v - 1], m_origdata->Column(v));
            if (m_printF)
                LOG << "SpatialVbInferenceTechnique::Fnoise == " << F << endl;

            //* MOVED HERE on Michael's advice -- 2007-11-23
            if (!m_locked_linear)
                m_lin_model[v - 1].ReCentre(m_fwd_post[v - 1].means);

            if (m_needF)
                F = noise->CalcFreeEnergy(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                    m_fwd_post[v - 1], m_fwd_prior[v - 1],
                    m_lin_model[v - 1], m_origdata->Column(v));
            if (m_printF)
                LOG << "SpatialVbInferenceTechnique::Flin == " << F << endl;
        }

        // next iteration:
        ++it;
    } while (!m_conv->Test(globalF));

    // Phew!

    // Interesting addition: calculate "coefficient resels" from Penny et al. 2005
    for (int k = 1; k <= m_num_params; k++)
    {
        ColumnVector gamma_vk(m_nvoxels); // might be handy
        ColumnVector gamma_vk_eo(
            m_nvoxels); // slightly different calculation (differs if using EO)
        gamma_vk_eo = -999;
        for (int v = 1; v <= m_nvoxels; v++)
        {
            gamma_vk(v) = 1 - m_fwd_post[v - 1].GetCovariance()(k, k) / m_fwd_prior[v - 1].GetCovariance()(k, k);
            if (m_fwd_post_no_prior.at(v - 1) != NULL)
            {
                gamma_vk_eo(v) = m_fwd_post[v - 1].GetCovariance()(k, k) / m_fwd_post_no_prior[v - 1]->GetCovariance()(k, k);
            }
        }
        LOG << "SpatialVariationalBayes::Coefficient resels per voxel for param "
            << k << ": " << gamma_vk.Sum() / m_nvoxels << " (vb) or "
            << gamma_vk_eo.Sum() / m_nvoxels << " (eo)" << endl;
    }

    for (int v = 1; v <= m_nvoxels; v++)
    {
        resultMVNs[v - 1] = new MVNDist(m_fwd_post[v - 1], m_noise_post[v - 1]->OutputAsMVN());

        if (m_save_without_prior)
        {
            resultMVNsWithoutPrior.at(v - 1) = new MVNDist(
                *m_fwd_post_no_prior[v - 1], m_noise_post[v - 1]->OutputAsMVN());
            // Should probably save the noiseWithoutPriors, but don't need that yet
            // (ever?)
        }
    }

    // resultFs are already stored as we go along.

    if (!m_needF)
    {
        for (int v = 1; v <= m_nvoxels; v++)
            assert(resultFs.at(v - 1) == 9999);
        // check we're not throwing away anything useful

        resultFs.clear();
        // clearing resultFs here should prevent an F image from being saved.
    }

    // Save Sinvs if possible
    if (m_save_spatial_priors)
    {
        // Copied from MVNDist::Save.  There are enough subtle differences
        // to justify duplicating the code here.

        // FIXME this will not work right now as the data is not per-voxel
        LOG << "SpatialVariationalBayes::Not saving spatial priors as not implemented" << endl;
#if 0
		Matrix vols;

		vols.ReSize(m_num_params, m_nvoxels * m_nvoxels);
		for (int k = 1; k <= m_num_params; k++)
		{
			assert(Sinvs.at(k-1).Nrows() == m_nvoxels);
			Matrix full = Sinvs.at(k - 1); // easier to visualize if in full form
			vols.Row(k) = full.AsColumn().t();
		}

		volume4D<float> output(m_nvoxels, m_nvoxels, 1, m_num_params);
		output.set_intent(NIFTI_INTENT_SYMMATRIX, 1, 1, 1);
		output.setdims(1, 1, 1, 1);
		output.setmatrix(vols);
		// no unThresholding needed
		LOG << vols.Nrows() << "," << vols.Ncols() << endl;

		save_volume4D(output, outputDir + "/finalSpatialPriors");
        m_origdata->SaveVoxelData("finalSpatialPriors", vols, NIFTI_INTENT_SYMMATRIX);
#endif
    }

    // Delete stuff (avoid memory leaks)
    for (int v = 1; v <= m_nvoxels; v++)
    {
        delete m_noise_post[v - 1];
        delete m_noise_prior[v - 1];
        delete m_fwd_post_no_prior.at(v - 1);
    }
}

// Binary search for data(index) == num
// Assumes data is sorted ascending!!
// Either returns an index such that data(index) == num
//   or -1 if num is not present in data.
inline int binarySearch(const ColumnVector &data, int num)
{
    int first = 1, last = data.Nrows();

    while (first <= last)
    {
        int test = (first + last) / 2;

        if (data(test) < num)
        {
            first = test + 1;
        }
        else if (data(test) > num)
        {
            last = test - 1;
        }
        else if (data(test) == num)
        {
            return test;
        }
        else
        {
            assert(false); // logic error!  data wasn't sorted?
        }
    }
    return -1;
}

// static int sign(double d1)
//{
//	return (double(0) < d1) - (d1 < double(0));
//}

/**
 * Check voxels are listed in order
 *
 * Order must be increasing in z value, or if same
 * increasing in y value, and if y and z are same
 * increasing in x value.
 */
bool IsCoordMatrixCorrectlyOrdered(const Matrix &voxelCoords)
{
    // Only 3D
    assert(voxelCoords.Nrows() == 3);

    // Voxels are stored one per column, each column is the x/y/z coords
    const int m_nvoxels = voxelCoords.Ncols();

    // Go through each voxel one at a time apart from last
    for (int v = 1; v <= m_nvoxels - 1; v++)
    {
        // Find difference between current coords and next
        ColumnVector diff = voxelCoords.Column(v + 1) - voxelCoords.Column(v);

        // Check order
        // +1 = +x, +10 = +y, +100 = +z, -99 = -z+x, etc.
        int d = sign(diff(1)) + 10 * sign(diff(2)) + 100 * sign(diff(3));
        if (d <= 0)
        {
            //			LOG << "Found mis-ordered voxels " << v << " and
            //" << v + 1 << ": d=" << d << endl;
            return false;
        }
    }
    return true;
}

/**
 * Calculate nearest and second-nearest neighbours for the voxels
 *
 * FIXME should use member variable for co-ords but test needs fixing in this
 * case
 */
void SpatialVariationalBayes::CalcNeighbours(const Matrix &voxelCoords)
{
    m_coords = &voxelCoords; // FIXME temp hack for testing
    const int nVoxels = voxelCoords.Ncols();
    if (nVoxels == 0)
        return;

    // Voxels must be ordered by increasing z, y and x values respectively in
    // order
    // of priority otherwise binary search for voxel by offset will not work
    if (!IsCoordMatrixCorrectlyOrdered(voxelCoords))
    {
        throw FabberInternalError(
            "Coordinate matrix must be in correct order to "
            "use adjacency-based priors.");
    }

    // Create a column vector with one entry per voxel.
    ColumnVector offsets(nVoxels);

    // Populate offsets with the offset into the
    // matrix of each voxel. We assume that co-ordinates
    // could be zero but not negative
    int xsize = m_coords->Row(1).Maximum() + 1;
    int ysize = m_coords->Row(2).Maximum() + 1;
    for (int v = 1; v <= nVoxels; v++)
    {
        int x = (*m_coords)(1, v);
        int y = (*m_coords)(2, v);
        int z = (*m_coords)(3, v);
        int offset = z * xsize * ysize + y * xsize + x;
        offsets(v) = offset;
    }

    // Delta is a list of offsets to find nearest
    // neighbours in x y and z direction (not diagonally)
    // Of course applying these offsets naively would not
    // always work, e.g. offset of -1 in the x direction
    // will not be a nearest neighbour for the first voxel
    // so need to check for this in subsequent code
    vector<int> delta;
    delta.push_back(1);              // next row
    delta.push_back(-1);             // prev row
    delta.push_back(xsize);          // next column
    delta.push_back(-xsize);         // prev column
    delta.push_back(xsize * ysize);  // next slice
    delta.push_back(-xsize * ysize); // prev slice

    // Don't look for neighbours in all dimensions.
    // For example if spatialDims=2, max_delta=3 so we
    // only look for neighbours in rows and columns
    //
    // However note we still need the full list of 3D deltas for later
    int max_delta = m_spatial_dims * 2 - 1;

    // Neighbours is a vector of vectors, so each voxel
    // will have an entry which is a vector of its neighbours
    m_neighbours.resize(nVoxels);

    // Go through each voxel. Note that offsets is indexed from 1 not 0
    // however the offsets themselves (potentially) start at 0.
    for (int vid = 1; vid <= nVoxels; vid++)
    {
        // Get the voxel offset into the matrix
        int pos = int(offsets(vid));

        // Now search for neighbours
        for (int n = 0; n <= max_delta; n++)
        {
            // is there a voxel at this neighbour position?
            // indexed from 1; id == -1 if not found.
            int id = binarySearch(offsets, pos + delta[n]);

            // No such voxel: continue
            if (id < 0)
                continue;

            // Check for wrap-around

            // Don't check for wrap around on final co-ord
            // PREVIOUSLY		if (delta.size() >= n + 2)
            // Changed (fixed)? because if spatialDims != 3 we still need
            // to check for wrap around in y-coordinate FIXME check
            if (n < 4)
            {
                bool ignore = false;
                if (delta[n] > 0)
                {
                    int test = delta[n + 2];
                    if (test > 0)
                        ignore = (pos % test) >= test - delta[n];
                }
                else
                {
                    int test = -delta[n + 2];
                    if (test > 0)
                        ignore = (pos % test) < -delta[n];
                }
                if (ignore)
                {
                    continue;
                }
            }
            // If we get this far, add it to the list
            m_neighbours.at(vid - 1).push_back(id);
            // cout << "voxel " << vid << " has neighbour " << id << endl;
        }
    }

    // Similar algorithm but looking for Neighbours-of-neighbours, excluding self,
    // but apparently including duplicates if there are two routes to get there
    // (diagonally connected) FIXME is this behaviour intentional?
    m_neighbours2.resize(nVoxels);

    for (int vid = 1; vid <= nVoxels; vid++)
    {
        // Go through the list of neighbours for each voxel.
        for (unsigned n1 = 0; n1 < m_neighbours.at(vid - 1).size(); n1++)
        {
            // n1id is the voxel index (not the offset) of the neighbour
            int n1id = m_neighbours[vid - 1].at(n1);
            int checkNofN = 0;
            // Go through each of it's neighbours. Add each, apart from original voxel
            for (unsigned n2 = 0; n2 < m_neighbours.at(n1id - 1).size(); n2++)
            {
                int n2id = m_neighbours[n1id - 1].at(n2);
                if (n2id != vid)
                {
                    m_neighbours2[vid - 1].push_back(n2id);
                }
                else
                    checkNofN++;
            }

            if (checkNofN != 1)
            {
                throw FabberInternalError(
                    "Each of this voxel's neighbours must have "
                    "this voxel as a neighbour");
            }
        }
    }
}

#if 0
/**
 * Convert a mask specifying voxels of interest into a list of xyz co-ordinates of
 * voxels of interest
 *
 * @param mask Input mask volume whose entries are 1 for voxels of interest, 0 otherwise
 * @param voxelCoords Return matrix in which every column is a the xyz vector or a
 *                    voxel's co-ordinates
 */
void ConvertMaskToVoxelCoordinates(const volume<float>& mask, Matrix& voxelCoords)
{
// Mask has previously been binarized to 0 or 1 to identify
// voxels of interest, so sum is count of these voxels
	ColumnVector preThresh((int) mask.sum());
	const int nVoxels = preThresh.Nrows();

// Populate preThresh with offset of voxel into matrix
// starting at 0 FIXME repeated from calculateNeigbours
	int offset(0);
	int count(1);
	for (int z = 0; z < mask.zsize(); z++)
	{
		for (int y = 0; y < mask.ysize(); y++)
		{
			for (int x = 0; x < mask.xsize(); x++)
			{
				if (mask(x, y, z) != 0)
				{
					preThresh(count++) = offset;
				}
				offset++;
			}
		}
	}

// Dimensions of mask matrix
	const int dims[3] =
	{	mask.xsize(), mask.ysize(), mask.zsize()};

// Basic sanity checks. Note that xdim/ydim/zdim are the physical
// dimensions of the voxels
	assert(mask.xsize()*mask.ysize()*mask.zsize() > 0);
	assert(mask.xdim()*mask.ydim()*mask.zdim() > 0);

	LOG_SAFE_ELSE_DISCARD("Calculating distance matrix, using voxel dimensions: "
			<< mask.xdim() << " by " << mask.ydim() << " by " << mask.zdim() << " mm\n" );

	ColumnVector positions[3];// indices
	positions[0].ReSize(nVoxels);
	positions[1].ReSize(nVoxels);
	positions[2].ReSize(nVoxels);

// Go through each voxel. Note that preThresh is a ColumnVector and indexes from 1 not 0
	for (int vox = 1; vox <= nVoxels; vox++)
	{
		int pos = (int) preThresh(vox) - 1; // preThresh appears to be 1-indexed. FIXME not sure, from above looks like 0 offset is possible for first voxel
		assert(pos>=0);

		positions[0](vox) = pos % dims[0];// X co-ordinate of offset from preThresh
		pos = pos / dims[0];
		positions[1](vox) = pos % dims[1];// Y co-ordinate of offset from preThresh
		pos = pos / dims[1];
		positions[2](vox) = pos % dims[2];// Z co-ordinate of offset from preThresh
		pos = pos / dims[2];
//LOG << vox << ' ' << (int)preThresh(vox)<< ' ' << dims[0] << ' ' << dims[1] << ' ' << dims[2] << ' ' << pos << endl;
		assert(pos == 0);// fails if preThresh >= product of dims[0,1,2]

// FIXME looks like old code below - remove?

//double pos = preThresh(vox);
//positions[2](vox) = floor(pos/dims[0]/dims[1]);
//pos -= positions[2](vox)*dims[0]*dims[1];
//positions[1](vox) = floor(pos/dims[0]);
//pos -= positions[1](vox)*dims[0];
//positions[0](vox) = pos;
//assert(positions[2](vox) < dims[2]);
//assert(positions[1](vox) < dims[1]);
//assert(positions[0](vox) < dims[0]);
		assert(preThresh(vox)-1 == positions[0](vox) + positions[1](vox)*dims[0] + positions[2](vox)*dims[0]*dims[1] );
	}

// Turn 3 column vectors into a matrix where each column is a voxel. note that
// need to transpose to do this
	voxelCoords = (positions[0] | positions[1] | positions[2]).t();
}
#endif

/**
 * Calculate a distance matrix
 *
 * FIXME voxelCoords should really be in MM, not indices; only really matters if
 * it's aniostropic or you're using the
 * smoothness values directly.
 *
 * @param voxelCoords List of voxel co-ordinates as a matrix: each column = 1
 * voxel
 * @param distanceMeasure How to measure distance: dist1 = Euclidian distance,
 * dist2 = squared Euclidian distance,
 *                        mdist = Manhattan distance (|dx| + |dy|)
 */
void CovarianceCache::CalcDistances(const NEWMAT::Matrix &voxelCoords,
    const string &distanceMeasure)
{
    // Create 3 column vectors, one to hold X co-ordinates, one Y and one Z
    ColumnVector positions[3];
    positions[0] = voxelCoords.Row(1).t();
    positions[1] = voxelCoords.Row(2).t();
    positions[2] = voxelCoords.Row(3).t();

    const int nVoxels = positions[0].Nrows();

    // dimSize is already included in voxelCoords
    // FIXME not obvious that it is, if not this should
    // be the dimensions of a voxel in mm
    const double dimSize[3] = { 1.0, 1.0, 1.0 };

    if (nVoxels > 7500)
    {
        LOG << "SpatialVariationalBayes::CalcDistances Over " << int(2.5 * nVoxels * nVoxels * 8 / 1e9)
            << " GB of memory will be used just to calculate "
            << "the distance matrix.  Hope you're not trying to invert this "
               "sucker!\n"
            << endl;
    }

    // 3 NxN symmetric matrices where N is the number of voxels
    // Each entry gives the absolute difference between
    // X/Y/Z co-ordinates respectively (FIXME in millimetres
    // supposedly but perhaps not given comments above)
    SymmetricMatrix relativePos[3];

    // Column vector with one entry per voxel, all set to 1 initially
    ColumnVector allOnes(nVoxels);
    allOnes = 1.0;

    // FIXME do not understand why this works!
    for (int dim = 0; dim < 3; dim++)
    {
        Matrix rel = dimSize[dim] * (positions[dim] * allOnes.t() - allOnes * positions[dim].t());
        assert(rel == -rel.t());
        assert(rel.Nrows() == nVoxels);
        // Down-convert to symmetric matrix (lower triangle so all positive??)
        relativePos[dim] << rel;
    }

    // Distances is an NxN symmetric matrix where N is number of voxels
    distances.ReSize(nVoxels);

    // FIXME code repetition here, distance measure should invoke a function
    // probably
    if (distanceMeasure == "dist1")
    {
        // absolute Euclidean distance
        LOG << "SpatialVariationalBayes::Using absolute Euclidean distance\n";
        for (int a = 1; a <= nVoxels; a++)
        {
            for (int b = 1; b <= a; b++)
            {
                distances(a, b) = sqrt(relativePos[0](a, b) * relativePos[0](a, b) + relativePos[1](a, b) * relativePos[1](a, b) + relativePos[2](a, b) * relativePos[2](a, b));
            }
        }
    }
    else if (distanceMeasure == "dist2")
    {
        // Euclidian distance squared
        LOG << "SpatialVariationalBayes::Using almost-squared (^1.99) Euclidean distance\n";
        for (int a = 1; a <= nVoxels; a++)
        {
            for (int b = 1; b <= a; b++)
            {
                distances(a, b) = pow(relativePos[0](a, b) * relativePos[0](a, b) + relativePos[1](a, b) * relativePos[1](a, b) + relativePos[2](a, b) * relativePos[2](a, b),
                    0.995);
            }
        }
    }
    else if (distanceMeasure == "mdist")
    {
        // Manhattan distance (bad?)
        LOG << "SpatialVariationalBayes::Using Manhattan distance\n";
        LOG << "SpatialVariationalBayes::WARNING: Seems to result in numerical problems down the line (not "
               "sure why)\n";
        for (int a = 1; a <= nVoxels; a++)
        {
            for (int b = 1; b <= a; b++)
            {
                distances(a, b) = fabs(relativePos[0](a, b)) + fabs(relativePos[1](a, b)) + fabs(relativePos[2](a, b));
            }
        }
    }
    else
    {
        throw InvalidOptionValue("distance-measure", distanceMeasure,
            "Unrecognized distance measure");
    }
}

// /*
class DerivFdRho : public GenericFunction1D
{
public:
    virtual double Calculate(double input) const;
    DerivFdRho(const CovarianceCache &cov,
        const DiagonalMatrix &cr,
        const ColumnVector &mdr,
        const double d)
        : m_covar(cov)
        , covRatio(cr)
        , meanDiffRatio(mdr)
        , delta(d)
    {
        return;
    }

    // don't need a PickFasterGuess

private:
    const CovarianceCache &m_covar;
    const DiagonalMatrix &covRatio;
    const ColumnVector &meanDiffRatio;
    const double delta;
};

double DerivFdRho::Calculate(const double rho) const
{
    const int Nvoxels = m_covar.GetDistances().Nrows();
    const SymmetricMatrix &Cinv = m_covar.GetCinv(delta);

    double out = 0;
    out += 0.5 * Nvoxels;
    out += -0.5 * (covRatio * exp(rho) * Cinv).Trace();
    out += -0.5 * (meanDiffRatio.t() * exp(rho) * Cinv * meanDiffRatio).AsScalar();

    /* Old version (pre Oct 10) -- actually gives almost-identical results (just
   the prior).

     double out = Nvoxels; // == Trace(S*Sinv)
     //  out -= exp(rho) * SP(covRatio, Cinv).Trace();
     out -= exp(rho) * (covRatioSupplemented * Cinv).Trace();
     out -= exp(rho) * (meanDiffRatio.t() * Cinv * meanDiffRatio).AsScalar();

     // Want a scale-free prior on exp(rho)

     //  LOG_ERR("Scale-free prior on rho... +exp(-2*rho)\n");
     out -= -1/exp(2*rho); // correct?

     */

    return out;
}
// */

// Evidence optimization
class DerivEdDelta : public GenericFunction1D
{
public:
    virtual double Calculate(double delta) const;
    DerivEdDelta(const CovarianceCache &c,
        const vector<MVNDist *> &fpwp,
        //	       const vector<SymmetricMatrix>& Si)
        const int kindex,
        const MVNDist *initFwdPrior,
        const bool r = false)
        : m_covar(c)
        , m_fwd_post_no_prior(fpwp)
        , k(kindex)
        , initialFwdPrior(initFwdPrior)
        , allowRhoToVary(r)
    // Sinvs(Si)
    {
        return;
    }

    //  virtual bool PickFasterGuess(double* guess, double lower, double upper,
    //  bool allowEndpoints = false) const
    virtual double OptimizeRho(double delta) const;

private:
    const CovarianceCache &m_covar; // stores C, Cinv, distance matrix?, etc.
    //  const SymmetricMatrix& XXtr; // = X*X'*precision -- replaces covRatio
    //  const ColumnVector& XYtr; // = X*Y'*precision -- replaces meanDiffRatio;
    const vector<MVNDist *> &m_fwd_post_no_prior;

    // To allow multiple parameters, will need storage for k (which param to
    // optimize)
    // and the other prior precision matrices (may be from Penny, Nonspatial,
    // etc.)

    const int k;
    const MVNDist *initialFwdPrior;

    // const vector<SymmetricMatrix>& Sinvs;

    bool allowRhoToVary;
};

double DerivEdDelta::OptimizeRho(double delta) const
{
    double rho;
    if (!allowRhoToVary)
        return 0.0;

    // This is just copy-pasted from ::Calculate.  There are more efficient
    // ways to do this!
    const SymmetricMatrix &dist = m_covar.GetDistances();
    const int Nvoxels = dist.Nrows();

    assert(initialFwdPrior->GetCovariance()(k, k) == 1); // unimplemented correction factor!

    DiagonalMatrix XXtr(Nvoxels);
    ColumnVector XYtr(Nvoxels);
    {
        assert(Nvoxels == (int)m_fwd_post_no_prior.size());
        for (int v = 1; v <= Nvoxels; v++)
        {
            XXtr(v, v) = m_fwd_post_no_prior[v - 1]->GetPrecisions()(k, k);
            XYtr(v) = XXtr(v, v) * (m_fwd_post_no_prior[v - 1]->means(k) - initialFwdPrior->means(k));
        }
        assert(XXtr.Nrows() == Nvoxels);
        assert(XYtr.Nrows() == Nvoxels);
    }

    SymmetricMatrix Sigma;
    {
        Sigma = (XXtr + m_covar.GetCinv(delta)).i();
    }

    const ColumnVector mu = Sigma * XYtr;

    rho = -log(1.0 / Nvoxels * ((Sigma + mu * mu.t()) * m_covar.GetCinv(delta)).Trace());

    LOG << "DerivEdDelta::OptimizeRho rho == " << rho << endl;

    return rho;
}

double DerivEdDelta::Calculate(double delta) const
{
    //  assert(delta >= 0.05); // Will be slow below this scale

    const SymmetricMatrix &dist = m_covar.GetDistances();
    const int Nvoxels = dist.Nrows();

    DiagonalMatrix XXtr(Nvoxels);
    ColumnVector XYtr(Nvoxels);

    {
        assert(Nvoxels == (int)m_fwd_post_no_prior.size());
        for (int v = 1; v <= Nvoxels; v++)
        {
            XXtr(v, v) = m_fwd_post_no_prior[v - 1]->GetPrecisions()(k, k) * initialFwdPrior->GetCovariance()(k, k);
            XYtr(v) = XXtr(v, v) * (m_fwd_post_no_prior[v - 1]->means(k) - initialFwdPrior->means(k)) * sqrt(initialFwdPrior->GetPrecisions()(k, k));
            WARN_ONCE("DerivEdDelta::Using the new XYtr correction (*sqrt(precision))");
        }
        assert(XXtr.Nrows() == Nvoxels);
        assert(XYtr.Nrows() == Nvoxels);
    }

    double out; // Assigned by the following statement:
    const SymmetricMatrix &CiCodistCi = m_covar.GetCiCodistCi(delta, &out);
    SymmetricMatrix Sigma;
    {
        Sigma = (XXtr + m_covar.GetCinv(delta)).i();
    }

    out -= (Sigma * CiCodistCi).Trace();

    // If the trace turns out to be slow, then
    // use identity: trace(a*b) == sum(sum(a.*b'))

    const ColumnVector mu = Sigma * XYtr;

    out -= (mu.t() * CiCodistCi * mu).AsScalar();
    out /= -4 * delta * delta; // = -1/2 * d(1/delta)/ddelta.  Note Sahani used
    // d(1/delta^2)/ddelta.

    if (0)
    {
        // WORKS! This is a very strong prior that attracts delta towards 5.
        //	LOG_ERR("Using a Ga(.0001,50000) prior on delta!\n");
        //	const double c = 50000, b = .0001;
        //	// DERIVATIVE OF -gammaln(c) + (c-1)*log(delta) - c*log(b) - delta/b;
        //	out += (c-1)/delta - 1/b;
        // Note: Matlab is gampdf(x, c, b) NOT b,c!

        // Uninformative prior should have b >>> 1.
        const double b = 1e40;
        const double c = 1 / b;
        // Mean = b*c
        // Variance = b^2 * c

        // Strong prior:
        //	const double m = 3.0/4.0; // correcting for 1/4 scaling of this ASL data
        //	const double s = 1.0/4.0;
        //
        //	const double b = (s*s)/m;
        //	const double c = (m*m)/(s*s);

        WARN_ONCE("DerivEdDelta::Using a Ga(" + stringify(b) + ", " + stringify(c) + " prior on delta!");
        // DERIVATIVE OF -gammaln(c) + (c-1)*log(delta) - c*log(b) - delta/b;
        out += (c - 1) / delta - 1 / b;
    }
    else
    {
        //	LOG_ERR("Not using any prior at all on delta\n");
        WARN_ONCE("DerivEdDelta::Not using any prior at all on delta");
    }

    return out;
}

// helper function
void BlockInverse(const SymmetricMatrix &in, SymmetricMatrix &out);

// Old Free energy optimization
class DerivFdDelta : public GenericFunction1D
{
public:
    virtual double Calculate(double input) const;
    DerivFdDelta(const CovarianceCache &cov,
        const DiagonalMatrix &cr,
        const ColumnVector &mdr,
        bool rv = true)
        : m_covar(cov)
        , covRatio(cr)
        , meanDiffRatio(mdr)
        , allowRhoToVary(rv)
    {
        return;
    }

    virtual bool PickFasterGuess(double *guess,
        double lower,
        double upper,
        bool allowEndpoints = false) const
    {
        return m_covar.GetCachedInRange(guess, lower, upper, allowEndpoints);
    }

    virtual double OptimizeRho(double delta) const
    {
        if (!allowRhoToVary)
            return 0.0;

        // Use the following approximation, without the prior.
        // The difference due to prior is typically very small..
        // >5: identical to 6 digits
        // 1.80468 -> 1.80471
        // -7 -> -4
        // -12 -> -6
        // -14 -> -7
        // -19 -> -8
        // In my data set, these correspond to delta = 1e12 or 1e15
        // For values with a rho (not dt), typically <1000 and highest
        // observed stop value was 700,000.

        const int Nvoxels = m_covar.GetDistances().Nrows();
        const SymmetricMatrix &Cinv = m_covar.GetCinv(delta);
        //      const double tmp = SP(covRatio, Cinv).Trace()
        const double tmp = (covRatio * Cinv).Trace() + (meanDiffRatio.t() * Cinv * meanDiffRatio).AsScalar();
        // Note: tmp can be negative if there's a numerical problem.
        // this means rho2 = NaN so it'll go on to the search method,
        // which should deal with this case reasonably well...

        //      LOG << Cinv.Trace()<<endl
        //	   << SP(Cinv, IdentityMatrix(Nvoxels)) << endl
        //<< covRatio
        // << meanDiffRatio.t()
        //	   << SP(covRatio, Cinv).Trace() << endl
        //	   << (meanDiffRatio.t() * Cinv * meanDiffRatio).AsScalar()
        //	   << "  tmp == " << tmp << endl;

        const double rho2 = -log(tmp / Nvoxels);
        LOG << "DerivFddDelta::OptimizeRho rho2 == " << rho2 << endl;

        double rho;
        if (true) //(rho2 > 2.0)
        {
            rho = rho2;
        }
        else
        {
            // Rho is small enough that the prior might actually
            // make a difference -- so calculate it the more-accurate-but-slow
            // way.

            //	  LOG_ERR("\n--- OPTIMIZING FOR RHO ---\n");
            DerivFdRho fcn2(m_covar, covRatio, meanDiffRatio, delta);
            fcn2.SetLogger(m_log);
            BisectionGuesstimator guesser;
            guesser.SetLogger(m_log);
            rho = DescendingZeroFinder(fcn2)
                      .InitialGuess(1)
                      .TolY(0.0001)
                      .RatioTolX(1.001)
                      .Verbosity(0)
                      .SetGuesstimator(&guesser)
                      .SearchMin(-70)
                      .SearchMax(70);

            // Observed range: -35 to +40
            // If the inversion causes numerical problems, then the above
            // we get tmp < 0, so rho2 = NaN, so we search, and typically end
            // up stuck at the upper limit.  700 seems to cause problems, but it
            // seems to recover if we only go as high as 70.

            //	  LOG_ERR("--- OPTIMIZED, RHO == " << rho << " ---\n\n");
            LOG << "DerivFddDelta::OptimizeRho rho == " << rho << endl;
        }

        return rho;
    }

private:
    const CovarianceCache &m_covar;
    const DiagonalMatrix &covRatio;
    // const SymmetricMatrix& covRatioSupplemented;
    const ColumnVector &meanDiffRatio;
    bool allowRhoToVary;
};

double DerivFdDelta::Calculate(const double delta) const
{
    const double rho = OptimizeRho(delta);
    // Returns rho = 0 if !allowRhoToVary.
    // Otherwise we return the value of delta optimized over all rhos!

    // Now that we've optimized over rho, get on with the calculation:

    //    LOG_ERR("        Calculating dF/ddelta at " << delta << endl);

    assert(delta >= 0.05);
    //    const SymmetricMatrix& dist = m_covar.GetDistances();
#ifndef NDEBUG
    const SymmetricMatrix &dist = m_covar.GetDistances();
    const int Nvoxels = dist.Nrows();
    assert(covRatio.Nrows() == Nvoxels);
    assert(meanDiffRatio.Nrows() == Nvoxels);
#endif

    // SymmetricMatrix C = m_covar.GetC(delta);
    // const SymmetricMatrix& Ci = m_covar.GetCinv(delta);
    // const Matrix& CiCodist = m_covar.GetCiCodist(delta);
    // double out = m_covar.GetCiCodist(delta).Trace();

    double out;
    const SymmetricMatrix &CiCodistCi = m_covar.GetCiCodistCi(delta, &out);
    // Above does: out = trace(CiCodist)

    //    LOG << "The uncacheable parts... " << flush;
    // LOG_ERR("values: " << out);

    // Correct???, but slower:
    //    out -= exp(rho) * (
    //		         (covRatio.i() - diag(old Ci) + Ci).i() * CiCodistCi
    //		      ).Trace();

    // Hopefully also correct (after iterations) but faster:
    // METHOD USED BEFORE 2008-03-13
    out -= exp(rho) * (covRatio * CiCodistCi).Trace();

    // If the trace turns out to be slow, then
    // use identity: trace(a*b) == sum(sum(a.*b'))

    out -= exp(rho) * (meanDiffRatio.t() * CiCodistCi * meanDiffRatio).AsScalar(); // REQUIRES MEANDIFFRATIO
    out /= -4 * delta * delta;
    //    LOG << "done." << endl;

    // Slower (equivalent) method
    // C = exp(-0.5*dist/delta_s);
    // Ci = inv(C);
    // Codist = C.*dist;
    // out2 = -0.25/delta_s^2 * ( ...
    //    trace(Ci * Codist) ...
    //    + trace(inv(diag(precRatio)) * Ci * Codist * Ci) ...
    //    + meanDiffRatio' * Ci * Codist * Ci * meanDiffRatio );
    //[out out2];
    // assert(diff(ans)/mean(ans) < 1e-10)

    //    LOG_ERR("        dF/ddelta at\t" << delta << " =\t" << out << endl);

    //    LOG_ERR(" - " << SP(covRatio, CiCodistCi).Trace()
    //    << " - " << (meanDiffRatio.t() * CiCodistCi * meanDiffRatio).AsScalar()
    //    << " all / " <<  -4*delta*delta
    //    << " - " << log(delta)
    //    <<" == " << out << endl);

    // Use a scale-free prior: F -= 1/delta
    //   so F' -= -1/delta^2

    //    out -= -1/delta/delta;
    //    LOG_ERR("Using NO PRIOR on delta\n");

    if (0)
    {
        // WORKS! This is a very strong prior that attracts delta towards 5.
        //	LOG_ERR("Using a Ga(.0001,50000) prior on delta!\n");
        //	const double c = 50000, b = .0001;
        //	// DERIVATIVE OF -gammaln(c) + (c-1)*log(delta) - c*log(b) - delta/b;
        //	out += (c-1)/delta - 1/b;

        WARN_ONCE("DerivFdDelta::Using a Ga(.1,50) prior on delta!");
        const double c = 50, b = .1;
        // DERIVATIVE OF -gammaln(c) + (c-1)*log(delta) - c*log(b) - delta/b;
        out += (c - 1) / delta - 1 / b;
    }
    else
    {
        WARN_ONCE("DerivFdDelta::Not using any prior at all on delta");
    }

    //    LOG_ERR("Using CORRECT scale-free prior on DerivFdDelta\n");

    return out;
}

double SpatialVariationalBayes::OptimizeEvidence(
    const vector<MVNDist *> &m_fwd_post_no_prior, // used for parameter k
    int k,
    const MVNDist *initialFwdPrior,
    double guess,
    bool allowRhoToVary,
    double *rhoOut) const
{
#ifndef NDEBUG
    const int Nparams = m_fwd_post_no_prior[0]->GetSize();
    assert(m_fwd_post_no_prior.at(0) != NULL);
    assert(Nparams >= 1);
    assert(k <= Nparams);
#endif

    DerivEdDelta fcn(m_covar, m_fwd_post_no_prior, k, initialFwdPrior,
        allowRhoToVary);
    fcn.SetLogger(m_log);

    LogBisectionGuesstimator guesser;
    guesser.SetLogger(m_log);

    double hardMin = 0.05; // Much below 0.2 inversion becomes painfully slow
    double hardMax = 1e3;  // Above 1e15, exp(-0.5*1/delta) == 1 (singular)

    double delta = DescendingZeroFinder(fcn)
                       .InitialGuess(guess)
                       .InitialScale(
                           guess * 0.009) // So hopefully it'll only take two guesses once settled
                       .ScaleGrowth(16)   // But it can escape pretty quickly!  7 guesses to
                       // get from 0.05 to 1000.
                       .SearchMin(hardMin)
                       .SearchMax(hardMax)
                       .RatioTolX(1.01)
                       .MaxEvaluations(2 + m_keep_param_covars)
                       .SetGuesstimator(&guesser);

    WARN_ONCE("SpatialVariationalBayes::OptimizeEvidence Hard limits on delta: [" + stringify(hardMin) + ", " + stringify(hardMax) + "]");

    if (rhoOut != NULL)
        *rhoOut = fcn.OptimizeRho(delta); // 0 if allowRhoToVary == false

    return delta;
}

double SpatialVariationalBayes::OptimizeSmoothingScale(
    const DiagonalMatrix &covRatio,
    const ColumnVector &meanDiffRatio,
    double guess,
    double *optimizedRho,
    bool allowRhoToVary,
    bool allowDeltaToVary) const
{
    DerivFdDelta fcn(m_covar, covRatio, meanDiffRatio, allowRhoToVary);
    fcn.SetLogger(m_log);

    LogBisectionGuesstimator guesser;
    guesser.SetLogger(m_log);

    if (m_brute_force_delta_search)
    {
        LOG << "SpatialVariationalBayes::BEGINNING BRUTE-FORCE DELTA SEARCH" << endl;
        LOG << "SpatialVariationalBayes::PARAMETERS:\ncovRatio = [" << covRatio << endl
            << "];\nmeanDiffRatio = [" << meanDiffRatio << "];\n";

        // BELOW: changed cutoff (was 1e16 originally)
        for (double dk = 0.001; dk < 1e4; dk *= (sqrt((2))))
        {
            LOG << "SpatialVariationalBayes::dk = " << dk << endl;
            LOG << "SpatialVariationalBayes::BRUTEFORCE=" << dk << "\t"
                << -0.5 * m_covar.GetC(dk).LogDeterminant().LogValue() << "\t"
                << -0.5 * (m_covar.GetCinv(dk) * covRatio).Trace() << "\t"
                << -0.5 * (meanDiffRatio.t() * m_covar.GetCinv(dk) * meanDiffRatio).AsScalar()
                << endl;
        }
        LOG << "SpatialVariationalBayes::END OF BRUTE-FORCE DELTA SEARCH" << endl;
    }

    double delta;
    if (allowDeltaToVary)
    {
        delta = DescendingZeroFinder(fcn)
                    .InitialGuess(guess)
                    .SearchMin(0.2) // Below this, inversion becomes painfully slow
                    //.SearchMin(0.01) // Below this, inversion becomes painfully slow
                    .SearchMax(1e15) // Above this, exp(-0.5*1/delta) == 1 (singular)
                    .RatioTolX(1.01)
                    .MaxEvaluations(2 + m_keep_param_covars)
                    .SetGuesstimator(&guesser);
        // LOG_ERR("HORRIBLE HACK: delta *= 0.95;\n");
        //	delta *= 0.95;
    }
    else
    {
        delta = guess;
    }

    if (allowDeltaToVary && optimizedRho != NULL)
        *optimizedRho = fcn.OptimizeRho(delta);

    return delta;
}

const ReturnMatrix CovarianceCache::GetC(double delta) const
{
    const int Nvoxels = distances.Nrows();

    SymmetricMatrix C(Nvoxels);
    if (delta == 0)
    {
        C = IdentityMatrix(Nvoxels);
    }
    else
    {
        for (int a = 1; a <= Nvoxels; a++)
            for (int b = 1; b <= a; b++)
                C(a, b) = exp(-0.5 * distances(a, b) / delta);
    }

    // NOTE: when distances = squared distance, prior is equivalent to white
    // noise smoothed with a Gaussian with sigma^2 = 2*delta (haven't actually
    // double-checked this yet).
    // BEWARE: delta is measured in millimeters!! (based on NIFTI file info).

    C.Release();
    return C;
}

bool CovarianceCache::GetCachedInRange(double *guess,
    double lower,
    double upper,
    bool allowEndpoints) const
{
    assert(guess != NULL);
    const double initialGuess = *guess;
    if (!(lower < initialGuess && initialGuess < upper))
    {
        LOG << "SpatialVariationalBayes::Uh-oh... lower = " << lower << ", initialGuess = " << initialGuess
            << ", upper = " << upper << endl;
    }
    assert(lower < initialGuess && initialGuess < upper);

    Cinv_cache_type::iterator it = Cinv_cache.lower_bound(lower);
    if (it == Cinv_cache.end())
        return false;
    if (it->first == lower && !allowEndpoints)
        ++it;
    if (it == Cinv_cache.end())
        return false;
    if (it->first > upper)
        return false;
    if (it->first == upper && !allowEndpoints)
        return false;

    // Success -- we have at least one fast guess!
    *guess = it->first;

    //  LOG << "Found a guess! " << lower << " < " << *guess << " < " << upper <<
    //  endl;

    // Can we find a better one?
    while (++it != Cinv_cache.end() && it->first <= upper)
    {
        if (it->first == upper && !allowEndpoints)
            break;

        //      if ( abs(it->first - initialGuess) < abs(*guess - initialGuess) )
        if (it->first < initialGuess || it->first - initialGuess < initialGuess - *guess)
            *guess = it->first;

        //      LOG << "Improved guess! " << lower << " < " << *guess << " < " <<
        //      upper << endl;
    }

    assert(lower < *guess && *guess < upper);

    return true;
}

const SymmetricMatrix &CovarianceCache::GetCinv(double delta) const
{
#ifdef NOCACHE
    WARN_ONCE("CovarianceCache::GetCinv Cache is disabled to avoid memory problems!");
    cinv = GetC(delta);
    if (cinv.Nrows() > 0)
    {
        // Appears to be memory problem with inverting a 0x0 matrix
        cinv = cinv.i();
    }

    return cinv;
#else
    if (Cinv_cache[delta].Nrows() == 0)
    {
        //      LOG << "[" << flush;
        //      LOG << "GetCinv cache miss... " << flush;
        //      const int Nvoxels = distances.Nrows();
        //      SymmetricMatrix C(Nvoxels);
        //      for (int a = 1; a <= Nvoxels; a++)
        //	for (int b = 1; b <= a; b++)
        //	  C(a,b) = exp(-0.5*distances(a,b)/delta);
        //      Cinv_cache[delta] = C.i();
        Cinv_cache[delta] = GetC(delta).i();
        //      LOG << "done." << endl;
        //      LOG << "]" << flush;
    }
    else
    {
        //      LOG << "GetCinv cache hit!\n";
    }

    return Cinv_cache[delta];
#endif
}

const SymmetricMatrix &CovarianceCache::GetCiCodistCi(
    double delta,
    double *CiCodistTrace) const
{
    if (CiCodistCi_cache[delta].first.Nrows() == 0)
    {
#ifdef NOCACHE
        CiCodistCi_cache.clear();
#endif
        //      LOG << "{" << flush;
        GetCinv(delta); // for sensible messages, make sure cache hits
        // LOG << "GetCiCodistCi cache miss... " << flush;
        Matrix CiCodist = GetCinv(delta) * SP(GetC(delta), distances);
        CiCodistCi_cache[delta].second = CiCodist.Trace();
        Matrix CiCodistCi_tmp = CiCodist * GetCinv(delta);
        CiCodistCi_cache[delta].first << CiCodistCi_tmp; // Force symmetric

        { // check something
            double maxAbsErr = (CiCodistCi_cache[delta].first - CiCodistCi_tmp)
                                   .MaximumAbsoluteValue();
            if (maxAbsErr > CiCodistCi_tmp.MaximumAbsoluteValue() * 1e-5)
            // If that test fails, you're probably in trouble.
            // Reducing it to e.g. 1e-5 (to make dist2 work)
            //   => non-finite alpha in iteration 2
            // Reduced it to 1e-5 to get mdist to work....
            //   => same result.  (oops, was mdist2)
            // Reduced it to 1e-5 to get mdist to work, again...
            //   =>
            {
                LOG << "CovarianceCache::GetCiCodistCi matrix not symmetric!\nError = "
                    << maxAbsErr << ", maxabsvalue = "
                    << CiCodistCi_tmp.MaximumAbsoluteValue() << endl;
                assert(false);
            }
        }
        //      LOG << "}" << flush;
    }

    if (CiCodistTrace != NULL)
        (*CiCodistTrace) = CiCodistCi_cache[delta].second;
    return CiCodistCi_cache[delta].first;
}
