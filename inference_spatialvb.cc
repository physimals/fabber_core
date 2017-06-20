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
#include <iomanip>

using MISCMATHS::sign;

static OptionSpec OPTIONS[] = {
    { "spatial-dims", OPT_INT, "Number of spatial dimensions", OPT_NONREQ, "3" },
    { "spatial-speed", OPT_STR, "Restrict speed of spatial smoothing", OPT_NONREQ,
        "-1" },
    { "distance-measure", OPT_STR, "", OPT_NONREQ, "dist1" },
    { "param-spatial-priors", OPT_STR,
        "Type of spatial priors for each parameter, as a sequence of characters. "
        "N=nonspatial, M=Markov random field, P=Penny, A=ARD",
        OPT_NONREQ, "S+" },
    { "update-spatial-prior-on-first-iteration", OPT_BOOL, "", OPT_NONREQ, "" },
    { "locked-linear-from-mvn", OPT_MVN, "MVN file containing fixed centres for linearization", OPT_NONREQ,
        "" },
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
    args.Set("default-prior-type", "N");
    VariationalBayesInferenceTechnique::Initialize(fwd_model, args);

    m_prior_types_str = GetPriorTypesStr(m_prior_types);
    m_spatial_dims = GetSpatialDims(args);
    m_spatial_speed = args.GetDoubleDefault("spatial-speed", -1);
    assert(m_spatial_speed > 1 || m_spatial_speed == -1);

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
        case 'N':
        case 'I':
        case 'A':
            break;

        case 'm':
        case 'M':
        case 'p':
        case 'P':
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
    m_update_first_iter = args.GetBool("update-spatial-prior-on-first-iteration");

    // Locked linearizations, if requested
    m_locked_linear_file = args.GetStringDefault("locked-linear-from-mvn", "");
    m_locked_linear = (m_locked_linear_file != "");
}

void SpatialVariationalBayes::SetupPerVoxelDists(FabberRunData &allData)
{
    // Initialized in voxel loop below (from file or default as required)
    m_noise_post.resize(m_nvoxels, NULL);
    m_noise_prior.resize(m_nvoxels, NULL);
    m_fwd_post.resize(m_nvoxels);

    // Re-centred in voxel loop below
    m_lin_model.resize(m_nvoxels, LinearizedFwdModel(m_model));

    // Static initialization for all voxels currently
    m_fwd_prior.resize(m_nvoxels, *initialFwdPrior);

    // Loaded from file if required, otherwise initialized during calculation
    resultMVNs.resize(m_nvoxels, NULL);

    // Initialized during calculation
    resultFs.resize(m_nvoxels, 9999); // 9999 is a garbage default value

    // Whether to fix the linearization centres (default: false)
    vector<MVNDist *> lockedLinearDists;
    Matrix lockedLinearCentres;
    const int nNoiseParams = initialNoisePrior->OutputAsMVN().GetSize();
    if (m_locked_linear)
    {
        LOG << "SpatialVariationalBayes::Loading fixed linearization centres from the MVN '"
            << m_locked_linear_file << "'\nNOTE: This does not check if the correct "
                                   "number of parameters is present!\n";
        MVNDist::Load(lockedLinearDists, m_locked_linear_file, allData, m_log);
        lockedLinearCentres.ReSize(m_num_params, m_nvoxels);
    }

    if (m_continueFromFile != "")
    {
        LOG << "SpatialVbInferenceTechnique::Continuing from file "
            << m_continueFromFile << endl;
        InitMVNFromFile(m_continueFromFile, allData, paramFilename);
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

        m_noise_prior[v - 1] = initialNoisePrior->Clone();
        noise->Precalculate(*m_noise_post[v - 1], *m_noise_prior[v - 1],
            m_origdata->Column(v));
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
        
        switch (m_shrinkage_type)
        {
        case 'p':
        case 'P':
        case 'm':
        case 'M':
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
                else // P
                    tmp1 += sigmak(v, v) * ((nn + tiny) * (nn + tiny) + nn);
            }

            ColumnVector Swk = tiny * wk;

            for (int v = 1; v <= m_nvoxels; v++)
            {
                for (vector<int>::iterator v2It = m_neighbours[v - 1].begin();
                     v2It != m_neighbours.at(v - 1).end(); ++v2It)
                {
                    Swk(v) += wk(v) - wk(*v2It);
                }
                if (m_shrinkage_type == 'p' || m_shrinkage_type == 'm')
                    Swk(v) += wk(v) * (m_spatial_dims * 2 - m_neighbours.at(v - 1).size());
            }
            double tmp2 = Swk.SumSquare(); //(Swk.t() * Swk).AsScalar();

            if (m_shrinkage_type == 'm' || m_shrinkage_type == 'M')
                tmp2 = DotProduct(Swk, wk);

            LOG << "SpatialVariationalBayes::UpdateAkmean k=" << k << ", tmp1=" << tmp1 << ", tmp2=" << tmp2 << endl;
	    
            gk(k, k) = 1 / (0.5 * tmp1 + 0.5 * tmp2 + 0.1); // prior q1 == 10 (1/q1 == 0.1)
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

void SpatialVariationalBayes::SetFwdPriorShrinkageType(int v, const NEWMAT::DiagonalMatrix &akmean)
{
    const double tiny = 0; // turns out to be no longer necessary.
    double weight8 = 0;    // weighted +8
    ColumnVector contrib8(m_num_params);
    contrib8 = 0.0;
    for (vector<int>::iterator nidIt = m_neighbours[v - 1].begin();
         nidIt != m_neighbours[v - 1].end(); ++nidIt)
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
        assert(nn <= m_spatial_dims * 2);
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

    if (m_shrinkage_type == 'p' || m_shrinkage_type == 'm')
    {
        //	Penny-style DirichletBC priors -- ignoring initialFwdPrior completely!
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

    if (m_shrinkage_type == 'm') {
        // Dirichlet BCs on MRF
        mTmp = contrib8 / (8 * m_spatial_dims * 2); 
    }               
    if (m_shrinkage_type == 'M') {
        mTmp = contrib8 / (8 * (nn + 1e-8));
    }

    // equivalent, when non-spatial priors are very weak: m_fwd_prior[v-1].means = mTmp;
    m_fwd_prior[v - 1].means = m_fwd_prior[v - 1].GetCovariance() * (spatialPrecisions * mTmp + initialFwdPrior->GetPrecisions() * initialFwdPrior->means);

    if (m_shrinkage_type == 'm' || m_shrinkage_type == 'M')
        m_fwd_prior[v - 1].means = m_fwd_prior[v - 1].GetCovariance() * spatialPrecisions * mTmp; // = mTmp;
}

double SpatialVariationalBayes::SetFwdPrior(int v, bool isFirstIteration)
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
        }
        else {
            assert(false);
        }
    }

    assert(initialFwdPrior->GetPrecisions().Nrows() == spatialPrecisions.Nrows());
    // Should check that earlier!  It's possible for basis=1 and priors=2x2
    // to slip through.  TODO

    // Should check that initialFwdPrior was already diagonal --
    // this will cause real problems if it isn't!!
    // (Safe way: SP of the covariance matrices -- that'd force
    // diagonality while preserving individual variance.)
    DiagonalMatrix finalPrecisions = spatialPrecisions;
    
    ColumnVector finalMeans = priorMeans - spatialPrecisions.i() * weightedMeans;

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

    return Fard;
}

// M = Markov random field - normally used
// P = Alternative to M (Penny prior?)
// N = non-spatial prior (model default)
// A = ARD prior
// I = image prior
//
// m/p are variations on M/P with different edge behaviour (Dirichlet BCs)
// 
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

    // Make the neighbours[] lists if required
    if (m_prior_types_str.find_first_of("mMpP") != string::npos)
    {
        CalcNeighbours(*m_coords);
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
    akmean = 1e-8;

    // FIXME can't calculate free energy with spatial VB yet 
    // This value never changes
    const double globalF = 1234.5678;

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
            if (m_shrinkage_type != '-')
            {
                SetFwdPriorShrinkageType(v, akmean);
            }
            double Fard = SetFwdPrior(v, (it == 0));

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
                m_origdata->Column(v), NULL, 0);

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


        // ITERATE OVER VOXELS
        for (int v = 1; v <= m_nvoxels; v++) {
            double &F = resultFs.at(v - 1);

            PassModelData(v);

            noise->UpdateNoise(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                m_fwd_post[v - 1], m_lin_model[v - 1],
                m_origdata->Column(v));

            if (m_needF) {
                F = noise->CalcFreeEnergy(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                    m_fwd_post[v - 1], m_fwd_prior[v - 1],
                    m_lin_model[v - 1], m_origdata->Column(v));
                if (m_printF)
                    LOG << "SpatialVbInferenceTechnique::Fnoise == " << F << endl;
            }

            // MOVED HERE on Michael's advice -- 2007-11-23
            if (!m_locked_linear)
                m_lin_model[v - 1].ReCentre(m_fwd_post[v - 1].means);

            if (m_needF) {
                F = noise->CalcFreeEnergy(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                    m_fwd_post[v - 1], m_fwd_prior[v - 1],
                    m_lin_model[v - 1], m_origdata->Column(v));
                if (m_printF)
                    LOG << "SpatialVbInferenceTechnique::Flin == " << F << endl;
            }
        }

        ++it;
    } while (!m_conv->Test(globalF));

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
        }
        LOG << "SpatialVariationalBayes::Coefficient resels per voxel for param "
            << k << ": " << gamma_vk.Sum() / m_nvoxels << " (vb) or "
            << gamma_vk_eo.Sum() / m_nvoxels << " (eo)" << endl;
    }

    for (int v = 1; v <= m_nvoxels; v++)
    {
        resultMVNs[v - 1] = new MVNDist(m_fwd_post[v - 1], m_noise_post[v - 1]->OutputAsMVN());
    }

    // resultFs are already stored as we go along.

    if (!m_needF)
    {
        // check we're not throwing away anything useful
        for (int v = 1; v <= m_nvoxels; v++)
            assert(resultFs.at(v - 1) == 9999);

        // clearing resultFs here should prevent an F image from being saved.
        resultFs.clear();
    }

    // Delete stuff (avoid memory leaks)
    for (int v = 1; v <= m_nvoxels; v++)
    {
        delete m_noise_post[v - 1];
        delete m_noise_prior[v - 1];
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

void SpatialVariationalBayes::CheckCoordMatrixCorrectlyOrdered(const Matrix &coords)
{
    // Only 3D
    assert(coords.Nrows() == 3);

    // Voxels are stored one per column, each column is the x/y/z coords
    const int m_nvoxels = coords.Ncols();

    // Go through each voxel one at a time apart from last
    for (int v = 1; v <= m_nvoxels - 1; v++)
    {
        // Find difference between current coords and next
        ColumnVector diff = coords.Column(v + 1) - coords.Column(v);

        // Check order
        // +1 = +x, +10 = +y, +100 = +z, -100 = -z+x, etc.
        int d = sign(diff(1)) + 10 * sign(diff(2)) + 100 * sign(diff(3));
        if (d <= 0)
        {
            LOG << "Found mis-ordered voxels " << v << " and " << v + 1 << ": d=" << d << endl;
            throw FabberInternalError("Coordinate matrix must be in correct order to use adjacency-based priors.");
        }
    }
}

/**
 * Calculate nearest and second-nearest neighbours for the voxels
 */
void SpatialVariationalBayes::CalcNeighbours(const Matrix &coords)
{
    const int nVoxels = coords.Ncols();
    if (nVoxels == 0)
        return;

    // Voxels must be ordered by increasing z, y and x values respectively
    // otherwise binary search for voxel by offset will not work
    CheckCoordMatrixCorrectlyOrdered(coords);

    // Create a column vector with one entry per voxel.
    ColumnVector offsets(nVoxels);

    // Populate offsets with the offset into the
    // matrix of each voxel. We assume that co-ordinates
    // could be zero but not negative
    int xsize = coords.Row(1).Maximum() + 1;
    int ysize = coords.Row(2).Maximum() + 1;
    for (int v = 1; v <= nVoxels; v++)
    {
        int x = coords(1, v);
        int y = coords(2, v);
        int z = coords(3, v);
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
