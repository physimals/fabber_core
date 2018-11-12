/*  inference_vb.cc - Variational Bayes with optional spatial smoothing

 Adrian Groves and Matthew Webster, FMRIB Image Analysis Group

 Copyright (C) 2007-2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_vb.h"

#include "convergence.h"
#include "easylog.h"
#include "priors.h"
#include "run_context.h"
#include "tools.h"
#include "version.h"

#include <miscmaths/miscmaths.h>
#include <newmatio.h>

#include <math.h>

static OptionSpec OPTIONS[] = {
    { "noise", OPT_STR, "Noise model to use (white or ar1)", OPT_REQ, "" },
    { "convergence", OPT_STR, "Name of method for detecting convergence", OPT_NONREQ, "maxits" },
    { "max-iterations", OPT_STR,
        "number of iterations of VB to use with the maxits convergence detector", OPT_NONREQ,
        "10" },
    { "min-fchange", OPT_STR,
        "When using the fchange convergence detector, the change in F to stop at", OPT_NONREQ,
        "10" },
    { "max-trials", OPT_STR, "When using the trial mode convergence detector, the maximum number "
                             "of trials after an initial reduction in F",
        OPT_NONREQ, "10" },
    { "print-free-energy", OPT_BOOL, "Output the free energy", OPT_NONREQ, "" },
    { "mcsteps", OPT_INT, "Number of motion correction steps", OPT_NONREQ, "0" },
    { "continue-from-mvn", OPT_MVN, "Continue previous run from output MVN files", OPT_NONREQ, "" },
    { "output-only", OPT_BOOL, "Skip model fitting, just output requested data based on supplied "
                               "MVN. Can only be used with continue-from-mvn",
        OPT_NONREQ, "" },
    { "noise-initial-prior", OPT_MATRIX, "MVN of initial noise prior", OPT_NONREQ, "" },
    { "noise-initial-posterior", OPT_MATRIX, "MVN of initial noise posterior", OPT_NONREQ, "" },
    { "noise-pattern", OPT_STR, "repeating pattern of noise variances for each point (e.g. 12 "
                                "gives odd and even data points different variances)",
        OPT_NONREQ, "1" },
    { "PSP_byname<n>", OPT_STR, "Name of model parameter to use image prior", OPT_NONREQ, "" },
    { "PSP_byname<n>_type", OPT_STR, "Type of image prior to use for parameter <n> - I=image prior",
        OPT_NONREQ, "" },
    { "PSP_byname<n>_image", OPT_IMAGE, "Image prior for parameter <n>", OPT_NONREQ, "" },
    { "PSP_byname<n>_prec", OPT_FLOAT, "Precision to apply to image prior for parameter <n>",
        OPT_NONREQ, "" },
    { "PSP_byname<n>_transform", OPT_STR, "Transform to apply to parameter <n>", OPT_NONREQ, "" },
    { "allow-bad-voxels", OPT_BOOL,
        "Continue if numerical error found in a voxel, rather than stopping", OPT_NONREQ, "" },
    { "ar1-cross-terms", OPT_STR, "For AR1 noise, type of cross-linking (dual, same or none)",
        OPT_NONREQ, "dual" },
    { "spatial-dims", OPT_INT, "Number of spatial dimensions", OPT_NONREQ, "3" },
    { "spatial-speed", OPT_STR, "Restrict speed of spatial smoothing", OPT_NONREQ, "-1" },
    { "distance-measure", OPT_STR, "", OPT_NONREQ, "dist1" },
    { "param-spatial-priors", OPT_STR,
        "Type of spatial priors for each parameter, as a sequence of characters. "
        "N=nonspatial, M=Markov random field, P=Penny, A=ARD",
        OPT_NONREQ, "N+" },
    { "update-spatial-prior-on-first-iteration", OPT_BOOL, "", OPT_NONREQ, "" },
    { "locked-linear-from-mvn", OPT_MVN, "MVN file containing fixed centres for linearization",
        OPT_NONREQ, "" },
    { "" },
};

void Vb::GetOptions(vector<OptionSpec> &opts) const
{
    InferenceTechnique::GetOptions(opts);

    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string Vb::GetDescription() const
{
    return "Variational Bayes inference technique";
}
string Vb::GetVersion() const
{
    return fabber_version();
}
InferenceTechnique *Vb::NewInstance()
{
    return new Vb();
}
void Vb::Initialize(FabberRunData &rundata)
{
    InferenceTechnique::Initialize(rundata);

    // Figure out if F needs to be calculated every iteration
    m_saveF = rundata.GetBool("save-free-energy");
    m_printF = rundata.GetBool("print-free-energy");

    // Motion correction related setup - by default no motion correction
    m_num_mcsteps = convertTo<int>(rundata.GetStringDefault("mcsteps", "0"));

    m_spatial_dims = rundata.GetIntDefault("spatial-dims", 3, 0, 3);
    if (m_spatial_dims == 1)
    {
        WARN_ONCE("spatial-dims=1 is weird... hope you're just testing!");
    }
    else if (m_spatial_dims == 2)
    {
        WARN_ONCE("spatial-dims=2 may not work the way you expect");
    }

    // Locked linearizations, if requested
    m_locked_linear = rundata.GetStringDefault("locked-linear-from-mvn", "") != "";
}

void Vb::IgnoreVoxel(int v)
{
    LOG << "Vb::IgnoreVoxel This voxel will be ignored in further updates" << endl;

    m_ctx->ignore_voxels.push_back(v);

    // Remove voxel from lists of neighbours of other voxels.
    // We identify affected voxels by looking in the neighbour
    // lists for the bad voxel, because any voxel which has
    // the bad voxel as a neighbour will be a neighbour of the
    // bad voxel
    vector<int> nn = m_ctx->neighbours[v - 1];
    for (vector<int>::iterator i = nn.begin(); i != nn.end(); ++i)
    {
        // Reference to list of neighbours of some other voxel which
        // has the bad voxel as a neighbour
        vector<int> &n2 = m_ctx->neighbours[*i - 1];

        n2.erase(std::remove(n2.begin(), n2.end(), v), n2.end());
    }

    // Same for next-nearest-neighbours
    nn = m_ctx->neighbours2[v - 1];
    for (vector<int>::iterator i = nn.begin(); i != nn.end(); ++i)
    {
        // Reference to list of neighbours of some other voxel which
        // has the bad voxel as a neighbour
        vector<int> &n2 = m_ctx->neighbours2[*i - 1];

        n2.erase(std::remove(n2.begin(), n2.end(), v), n2.end());
    }
}

/**
 * Calculate free energy. Note that this is currently unused in spatial VB
 */
double Vb::CalculateF(int v, string label, double Fprior)
{
    double F = 1234.5678;
    if (m_needF)
    {
        F = m_ctx->m_noise->CalcFreeEnergy(*m_ctx->noise_post[v - 1], *m_ctx->noise_prior[v - 1],
            m_ctx->fwd_post[v - 1], m_ctx->fwd_prior[v - 1], m_ctx->m_lin_model[v - 1],
            m_ctx->m_origdata->Column(v));
        F += Fprior;
        m_ctx->resultFs[v - 1] = F;
        if (m_printF)
        {
            LOG << "Vb::F" << label << " = " << F << endl;
        }
    }
    return F;
}

void Vb::DebugVoxel(int v, const string &where)
{
    LOG << where << " - voxel " << v << " of " << m_ctx->nvoxels << endl;
    LOG << "Prior means: " << endl << m_ctx->fwd_prior[v - 1].means.t();
    LOG << "Prior precisions: " << endl << m_ctx->fwd_prior[v - 1].GetPrecisions();
    LOG << "Posterior means: " << endl << m_ctx->fwd_post[v - 1].means.t();
    LOG << "Noise prior means: " << endl << m_ctx->noise_prior[v - 1]->OutputAsMVN().means.t();
    LOG << "Noise prior precisions: " << endl
        << m_ctx->noise_prior[v - 1]->OutputAsMVN().GetPrecisions();
    LOG << "Centre: " << endl << m_ctx->m_lin_model[v - 1].Centre();
    LOG << "Offset: " << endl << m_ctx->m_lin_model[v - 1].Offset();
    LOG << "Jacobian: " << endl << m_ctx->m_lin_model[v - 1].Jacobian() << endl;
}

bool Vb::IsSpatial(FabberRunData &rundata) const
{
    if (rundata.GetString("method") == "spatialvb")
    {
        return true;
    }
    else
    {
        // Really clunky way to detect if any spatial priors have been specified
        vector<Parameter> params;
        m_ctx->m_model->GetParameters(rundata, params);
        for (vector<Parameter>::iterator iter = params.begin(); iter != params.end(); iter++)
        {
            switch (iter->prior_type)
            {
            case PRIOR_SPATIAL_M:
            case PRIOR_SPATIAL_m:
            case PRIOR_SPATIAL_P:
            case PRIOR_SPATIAL_p:
                return true;
            }
        }
    }
    return false;
}

void Vb::DoCalculations(FabberRunData &rundata)
{
    // extract data (and the coords) from rundata for the (first) VB run
    // Rows are volumes
    // Columns are (time) series
    // num Rows is size of (time) series
    // num Cols is size of volumes
    m_ctx = new RunContext();
    m_ctx->Initialize(rundata);
    
    // pass in some (dummy) data/coords here just in case the model relies upon it
    // use the first voxel values as our dummies FIXME this shouldn't really be
    // necessary, need to find way for model to know about the data beforehand.
    if (m_ctx->nvoxels > 0)
        m_ctx->PassModelData(1);

    if (rundata.GetBool("output-only"))
    {
        // Do no calculations - now we have set resultMVNs we can finish
        LOG << "Vb::DoCalculations output-only set - not performing any calculations" << endl;
    }
    else if (IsSpatial(rundata))
    {
        DoCalculationsSpatial(rundata);
    }
    else
    {
        DoCalculationsVoxelwise(rundata);
    }

    if (!m_needF)
    {
        // clearing resultFs here should prevent an F image from being saved.
        m_ctx->resultFs.clear();
    }
}

void Vb::DoCalculationsVoxelwise(FabberRunData &rundata)
{
    vector<Parameter> params;
    m_ctx->m_model->GetParameters(rundata, params);
    vector<Prior *> priors = PriorFactory(rundata).CreatePriors(params);

    // Loop over voxels
    for (int v = 1; v <= m_ctx->nvoxels; v++)
    {
        m_ctx->PassModelData(v);

        m_ctx->v = v;
        m_ctx->it = 0;

        // Save our model parameters in case we need to revert later.
        // Note need to save prior in case ARD is being used
        NoiseParams *const noisePosteriorSave = m_ctx->noise_post[v - 1]->Clone();
        MVNDist fwdPosteriorSave(m_ctx->fwd_post[v - 1]);
        MVNDist fwdPriorSave(m_ctx->fwd_prior[v - 1]);

        // Give an indication of the progress through the voxels;
        rundata.Progress(v, m_ctx->nvoxels);
        double F = 1234.5678;

        try
        {
            m_ctx->m_lin_model[v - 1].ReCentre(m_ctx->fwd_post[v - 1].means);
            m_ctx->m_conv[v - 1]->Reset();

            // START the VB updates and run through the relevant iterations (according to the
            // convergence testing)
            do
            {
                double Fprior = 0;

                if (m_ctx->m_conv[v - 1]->NeedRevert()) // revert to previous solution if the convergence
                                                 // detector calls for it
                {
                    *m_ctx->noise_post[v - 1] = *noisePosteriorSave;
                    m_ctx->fwd_post[v - 1] = fwdPosteriorSave;
                    m_ctx->fwd_prior[v - 1] = fwdPriorSave;
                    m_ctx->m_lin_model[v - 1].ReCentre(m_ctx->fwd_post[v - 1].means);
                    if (m_debug)
                        DebugVoxel(v, "Reverted");
                }

                // Save old values if called for
                if (m_ctx->m_conv[v - 1]->NeedSave())
                {
                    *noisePosteriorSave = *m_ctx->noise_post[v - 1]; // copy values, not pointer!
                    fwdPosteriorSave = m_ctx->fwd_post[v - 1];
                    fwdPriorSave = m_ctx->fwd_prior[v - 1];
                }

                for (int k = 0; k < m_ctx->m_num_params; k++)
                {
                    Fprior += priors[k]->ApplyToMVN(&m_ctx->fwd_prior[v - 1], *m_ctx);
                }

                if (m_debug)
                    DebugVoxel(v, "Applied priors");

                F = CalculateF(v, "before", Fprior);

                m_ctx->m_noise->UpdateTheta(*m_ctx->noise_post[v - 1], m_ctx->fwd_post[v - 1],
                    m_ctx->fwd_prior[v - 1], m_ctx->m_lin_model[v - 1], m_ctx->m_origdata->Column(v), NULL,
                    m_ctx->m_conv[v - 1]->LMalpha());

                if (m_debug)
                    DebugVoxel(v, "Updated params");

                F = CalculateF(v, "theta", Fprior);

                m_ctx->m_noise->UpdateNoise(*m_ctx->noise_post[v - 1], *m_ctx->noise_prior[v - 1],
                    m_ctx->fwd_post[v - 1], m_ctx->m_lin_model[v - 1], m_ctx->m_origdata->Column(v));

                if (m_debug)
                    DebugVoxel(v, "Updated noise");

                F = CalculateF(v, "phi", Fprior);

                // Linearization update
                // Update the linear model before doing Free energy calculation
                // (and ready for next round of theta and phi updates)
                m_ctx->m_lin_model[v - 1].ReCentre(m_ctx->fwd_post[v - 1].means);

                if (m_debug)
                    DebugVoxel(v, "Re-centered");

                F = CalculateF(v, "lin", Fprior);

                ++m_ctx->it;
            } while (!m_ctx->m_conv[v - 1]->Test(F));

            if (m_debug)
                LOG << "Converged after " << m_ctx->it << " iterations" << endl;

            // Revert to old values at last stage if required
            if (m_ctx->m_conv[v - 1]->NeedRevert())
            {
                *m_ctx->noise_post[v - 1] = *noisePosteriorSave;
                m_ctx->fwd_post[v - 1] = fwdPosteriorSave;
                m_ctx->fwd_prior[v - 1] = fwdPriorSave;
                m_ctx->m_lin_model[v - 1].ReCentre(m_ctx->fwd_post[v - 1].means);
            }

            delete noisePosteriorSave;
        }
        catch (FabberInternalError &e)
        {
            LOG << "Vb::Internal error for voxel " << v << " at " << m_ctx->m_coords->Column(v).t()
                << " : " << e.what() << endl;

            if (m_halt_bad_voxel)
                throw;
        }
        catch (NEWMAT::Exception &e)
        {
            LOG << "Vb::NEWMAT exception for voxel " << v << " at " << m_ctx->m_coords->Column(v).t()
                << " : " << e.what() << endl;

            if (m_halt_bad_voxel)
                throw;
        }

        // now write the results to resultMVNs
        try
        {
            m_ctx->resultMVNs.at(v - 1)
                = new MVNDist(m_ctx->fwd_post[v - 1], m_ctx->noise_post[v - 1]->OutputAsMVN());
            if (m_needF)
                m_ctx->resultFs.at(v - 1) = F;
        }
        catch (...)
        {
            // Even that can fail, due to results being singular
            LOG << "Vb::Can't give any sensible answer for this voxel; outputting zero +- "
                   "identity\n";
            MVNDist *tmp = new MVNDist(m_log);
            tmp->SetSize(m_ctx->fwd_post[v - 1].means.Nrows()
                + m_ctx->noise_post[v - 1]->OutputAsMVN().means.Nrows());
            tmp->SetCovariance(IdentityMatrix(tmp->means.Nrows()));
            m_ctx->resultMVNs.at(v - 1) = tmp;
            if (m_needF)
                m_ctx->resultFs.at(v - 1) = F;
        }
    }
    for (unsigned int i = 0; i < priors.size(); i++)
    {
        delete priors[i];
    }
}

void Vb::DoCalculationsSpatial(FabberRunData &rundata)
{
    // Pass in some (dummy) data/coords here just in case the model relies upon it
    // use the first voxel values as our dummies FIXME this shouldn't really be
    // necessary, need to find way for model to know about the data beforehand.
    if (m_ctx->nvoxels > 0)
        m_ctx->PassModelData(1);

    // Make the neighbours[] lists if required
    // if (m_prior_types_str.find_first_of("mMpP") != string::npos)
    if (true) // FIXME
    {
        m_ctx->CalcNeighbours(*m_ctx->m_coords, m_spatial_dims);
    }

    vector<Parameter> params;
    m_ctx->m_model->GetParameters(rundata, params);
    vector<Prior *> priors = PriorFactory(rundata).CreatePriors(params);

    // Spatial loop currently uses a global convergence detector FIXME
    // needs to change
    CountingConvergenceDetector conv;
    conv.Initialize(rundata);
    double Fglobal = 1234.5678;
    int maxits = convertTo<int>(rundata.GetStringDefault("max-iterations", "10"));

    // MAIN ITERATION LOOP
    do
    {
        LOG << endl << "*** Spatial iteration *** " << (m_ctx->it + 1) << endl;

        // Give an indication of the progress through the voxels;
        rundata.Progress(m_ctx->it, maxits);
        double Fprior = 0;

        // ITERATE OVER VOXELS
        for (int v = 1; v <= m_ctx->nvoxels; v++)
        {
            m_ctx->v = v;

            m_ctx->PassModelData(v);

            // The steps below are essentially the same as regular VB, although
            // the code looks different as the per-voxel dists are set up at the
            // start rather than as we go
            try
            {
                Fprior = 0;

                // Apply prior updates for spatial or ARD priors
                for (int k = 0; k < m_ctx->m_num_params; k++)
                {
                    Fprior += priors[k]->ApplyToMVN(&m_ctx->fwd_prior[v - 1], *m_ctx);
                }
                if (m_debug)
                    DebugVoxel(v, "Priors set");

                // Ignore voxels where numerical issues have occurred
                if (std::find(m_ctx->ignore_voxels.begin(), m_ctx->ignore_voxels.end(), v)
                    != m_ctx->ignore_voxels.end())
                {
                    LOG << "Ignoring voxel " << v << endl;
                    continue;
                }

                CalculateF(v, "before", Fprior);

                m_ctx->m_noise->UpdateTheta(*m_ctx->noise_post[v - 1], m_ctx->fwd_post[v - 1],
                    m_ctx->fwd_prior[v - 1], m_ctx->m_lin_model[v - 1], m_ctx->m_origdata->Column(v), NULL, 0);
                if (m_debug)
                    DebugVoxel(v, "Theta updated");

                CalculateF(v, "theta", Fprior);
            }
            catch (FabberInternalError &e)
            {
                LOG << "Vb::Internal error for voxel " << v << " at " << m_ctx->m_coords->Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
            catch (NEWMAT::Exception &e)
            {
                LOG << "Vb::NEWMAT exception for voxel " << v << " at " << m_ctx->m_coords->Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
        }

        Fglobal = 0;
        for (int v = 1; v <= m_ctx->nvoxels; v++)
        {
            try {
                // Ignore voxels where numerical issues have occurred
                if (std::find(m_ctx->ignore_voxels.begin(), m_ctx->ignore_voxels.end(), v)
                    != m_ctx->ignore_voxels.end())
                {
                    LOG << "Ignoring voxel " << v << endl;
                    continue;
                }

                m_ctx->PassModelData(v);

                m_ctx->m_noise->UpdateNoise(*m_ctx->noise_post[v - 1], *m_ctx->noise_prior[v - 1],
                    m_ctx->fwd_post[v - 1], m_ctx->m_lin_model[v - 1], m_ctx->m_origdata->Column(v));
                if (m_debug)
                    DebugVoxel(v, "Noise updated");

                CalculateF(v, "noise", Fprior);

                if (!m_locked_linear)
                {
                    m_ctx->m_lin_model[v - 1].ReCentre(m_ctx->fwd_post[v - 1].means);
                    if (m_debug)
                        DebugVoxel(v, "Re-centre");
                }
                
                Fglobal += CalculateF(v, "lin", Fprior);
            }
            catch (FabberInternalError &e)
            {
                LOG << "Vb::Internal error for voxel " << v << " at " << m_ctx->m_coords->Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
            catch (NEWMAT::Exception &e)
            {
                LOG << "Vb::NEWMAT exception for voxel " << v << " at " << m_ctx->m_coords->Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
        }

        ++m_ctx->it;
    } while (!conv.Test(Fglobal));

    // Interesting addition: calculate "coefficient resels" from Penny et al. 2005
    for (int k = 1; k <= m_ctx->m_num_params; k++)
    {
        ColumnVector gamma_vk(m_ctx->nvoxels);
        for (int v = 1; v <= m_ctx->nvoxels; v++)
        {
            gamma_vk(v) = 1
                - m_ctx->fwd_post[v - 1].GetCovariance()(k, k)
                    / m_ctx->fwd_prior[v - 1].GetCovariance()(k, k);
        }
        LOG << "Vb::Coefficient resels per voxel for param " << k << ": "
            << gamma_vk.Sum() / m_ctx->nvoxels << endl;
    }

    for (int v = 1; v <= m_ctx->nvoxels; v++)
    {
        m_ctx->resultMVNs[v - 1]
            = new MVNDist(m_ctx->fwd_post[v - 1], m_ctx->noise_post[v - 1]->OutputAsMVN());
    }
    for (unsigned int i = 0; i < priors.size(); i++)
    {
        delete priors[i];
    }
}

void Vb::SaveResults(FabberRunData &rundata) const
{
    InferenceTechnique::SaveResults(rundata);

    LOG << "Vb::Preparing to save results..." << endl;
    int nVoxels = m_ctx->resultMVNs.size();

    if (rundata.GetBool("save-noise-mean") | rundata.GetBool("save-noise-std"))
    {
        if (m_ctx->m_noise_params > 0)
        {
            LOG << "Vb::Writing noise" << endl;
            Matrix noiseMean, noiseStd;
            noiseMean.ReSize(m_ctx->m_noise_params, nVoxels);
            noiseStd.ReSize(m_ctx->m_noise_params, nVoxels);
            for (int vox = 1; vox <= nVoxels; vox++)
            {
                for (int i = 1; i <= m_ctx->m_noise_params; i++)
                {
                    noiseStd(i, vox) = sqrt(
                        m_ctx->resultMVNs[vox - 1]->GetCovariance()(i + m_ctx->m_num_params, i + m_ctx->m_num_params));
                    noiseMean(i, vox) = m_ctx->resultMVNs[vox - 1]->means(i + m_ctx->m_num_params);
                }
            }
            // FIXME was this being saved before? Should it be?
            if (rundata.GetBool("save-noise-mean"))
                rundata.SaveVoxelData("noise_means", noiseMean);
            if (rundata.GetBool("save-noise-std"))
                rundata.SaveVoxelData("noise_stdevs", noiseStd);
        }
    }

    // Save the Free Energy estimates
    if (m_saveF && !m_ctx->resultFs.empty())
    {
        LOG << "Vb::Writing free energy" << endl;
        assert((int)m_ctx->resultFs.size() == nVoxels);
        Matrix freeEnergy;
        freeEnergy.ReSize(1, nVoxels);
        for (int vox = 1; vox <= nVoxels; vox++)
        {
            freeEnergy(1, vox) = m_ctx->resultFs.at(vox - 1);
        }
        rundata.SaveVoxelData("freeEnergy", freeEnergy);
    }
    else
    {
        LOG << "Vb::Free energy wasn't recorded, so no freeEnergy data saved" << endl;
    }
    LOG << "Vb::Done writing results." << endl;
}
