/*  inference_vb.cc - VB inference technique class declarations

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_vb.h"

#include "convergence.h"
#include "motioncorr.h"

#include <newmat.h>

#include <algorithm>
#include <math.h>

using namespace NEWMAT;

std::ostream &operator<<(std::ostream &out, const PriorType &value)
{
    return out << "PriorType: Parameter " << value.m_param_name << " type: " << value.m_type << " filename: "
               << value.m_filename << " precision: " << value.m_prec << endl;
}

PriorType::PriorType()
{
}

PriorType::PriorType(unsigned int idx, vector<string> param_names, FabberRunData &data)
    : m_param_name(param_names[idx])
    , m_idx(idx)
    , m_prec(-1)
{
    // Complexity below is due to there being two ways of specifying
    // priors. One is using the param-spatial-priors option which is
    // a sequence of chars in model parameter order, one for each
    // parameter. A + character means 'use the previous value for all
    // remaining parameters'. An 'I' means an image prior and
    // the filename is specified separately using an image-prior<n> option
    string types = GetTypesString(data, param_names.size());
    m_type = data.GetStringDefault("default-prior-type", "-")[0];
    char current_type = '-';
    for (size_t i = 0; i < types.size(); i++)
    {
        if (i == idx)
        {
            if (types[i] == '+')
                m_type = current_type;
            else
                m_type = types[i];
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
    if (m_type == 'I')
    {
        m_filename = "image-prior" + stringify(idx + 1);
    }

    // Here is the second way of specifying priors which will override the above.
    // PSP_byname<n>=parameter name
    // PSP_byname<n>_type=character indicating type
    // PSP_byname<n>_image=image prior filename
    // PSP_byname<n>_prec=image prior precision
    int current = 0;
    while (true)
    {
        current++;
        string name = data.GetStringDefault("PSP_byname" + stringify(current), "stop!");
        if (name == "stop!")
            break;

        if (name == param_names[idx])
        {
            m_type = convertTo<char>(data.GetStringDefault("PSP_byname" + stringify(current) + "_type", stringify(m_type)));
            if (m_type == 'I')
            {
                m_filename = "PSP_byname" + stringify(current) + "_image";
            }
            m_prec = data.GetDoubleDefault("PSP_byname" + stringify(current) + "_prec", -1);
        }
    }

    if (m_type == 'I')
    {
        //		LOG << "PriorType::Reading Image prior (" << m_param_name << "): " << m_filename << endl;
        m_image = data.GetVoxelData(m_filename).AsRow();
    }
}

string PriorType::GetTypesString(FabberRunData &rundata, unsigned int num_params)
{
    string priors_str = rundata.GetStringDefault("param-spatial-priors", "");
    if (priors_str.size() > num_params)
    {
        // Have been given too many prior types
        //WARN_ONCE("--param-spatial-priors=" + priors_str + ", but there are only " + stringify(m_num_params) + " parameters!\n");
    }
    if (priors_str.size() < num_params)
    {
        // Expand '+' char, if present, to give correct number of parameters
        // If not, don't worry will just use default prior type for the missing
        //cerr << "Expanding param-spatial-priors from " << priors_str << endl;
        int deficit = num_params - priors_str.size();
        size_t plus_pos = priors_str.find("+");
        if (plus_pos != std::string::npos)
        {
            priors_str.insert(plus_pos, deficit, '+');
        }
        //cerr << "... to " << priors_str << endl;
    }

    return priors_str;
}

void PriorType::SetPrior(MVNDist *prior, int voxel)
{
    if (m_type == 'I')
    {
        prior->means(m_idx + 1) = m_image(voxel);
    }
    if (m_prec > 0)
    {
        SymmetricMatrix prec = prior->GetPrecisions();
        prec(m_idx + 1, m_idx + 1) = m_prec;
        prior->SetPrecisions(prec);
    }
}

static OptionSpec OPTIONS[] = {
    { "noise", OPT_STR, "Noise model to use (white or ar1)", OPT_REQ, "" },
    { "convergence", OPT_STR, "Name of method for detecting convergence", OPT_NONREQ, "maxits" },
    { "max-iterations", OPT_STR, "number of iterations of VB to use with the maxits convergence detector",
        OPT_NONREQ, "10" },
    { "min-fchange", OPT_STR, "When using the fchange convergence detector, the change in F to stop at", OPT_NONREQ,
        "10" },
    { "max-trials", OPT_STR,
        "When using the trial mode convergence detector, the maximum number of trials after an initial reduction in F",
        OPT_NONREQ, "10" },
    { "print-free-energy", OPT_BOOL, "Output the free energy", OPT_NONREQ, "" },
    { "mcsteps", OPT_INT, "Number of motion correction steps", OPT_NONREQ, "0" },
    { "continue-from-mvn", OPT_MVN, "Continue previous run from output MVN files", OPT_NONREQ, "" },
    { "output-only", OPT_BOOL,
        "Skip model fitting, just output requested data based on supplied MVN. Can only be used with continue-from-mvn",
        OPT_NONREQ, "" },
    { "fwd-initial-prior", OPT_MATRIX, "MVN of initial model prior. Important for spatial VB using D prior", OPT_NONREQ, "" },
    { "fwd-initial-posterior", OPT_MATRIX, "MVN of initial model posterior", OPT_NONREQ, "" },
    { "noise-initial-prior", OPT_MATRIX, "MVN of initial noise prior", OPT_NONREQ, "" },
    { "noise-initial-posterior", OPT_MATRIX, "MVN of initial noise posterior", OPT_NONREQ, "" },
    { "noise-pattern", OPT_STR,
        "repeating pattern of noise variances for each point (e.g. 12 gives odd and even data points different variances)",
        OPT_NONREQ, "1" },
    { "PSP_byname<n>", OPT_STR, "Name of model parameter to use image prior", OPT_NONREQ, "" },
    { "PSP_byname<n>_type", OPT_STR, "Type of image prior to use for parameter <n> - I=image prior",
        OPT_NONREQ, "" },
    { "PSP_byname<n>_image", OPT_IMAGE, "Image prior for parameter <n>", OPT_NONREQ, "" },
    { "PSP_byname<n>_prec", OPT_FLOAT, "Precision to apply to image prior for parameter <n>", OPT_NONREQ, "" },
    { "locked-linear-from-mvn", OPT_MVN, "MVN file containing fixed centres for linearization", OPT_NONREQ,
        "" },
    { "allow-bad-voxels", OPT_BOOL, "Continue if numerical error found in a voxel, rather than stopping",
        OPT_NONREQ, "" },
    { "ar1-cross-terms", OPT_STR, "For AR1 noise, type of cross-linking (dual, same or none)", OPT_NONREQ,
        "dual" },
    { "" },
};

InferenceTechnique *VariationalBayesInferenceTechnique::NewInstance()
{
    return new VariationalBayesInferenceTechnique();
}

void VariationalBayesInferenceTechnique::GetOptions(vector<OptionSpec> &opts) const
{
    InferenceTechnique::GetOptions(opts);
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

std::string VariationalBayesInferenceTechnique::GetDescription() const
{
    return "Variational Bayes inference technique. See Chappell et al IEEE Trans Sig Proc 57:1 (2009)";
}

string VariationalBayesInferenceTechnique::GetVersion() const
{
    return "1.0"; // FIXME
}

void VariationalBayesInferenceTechnique::InitializeMVNFromParam(FabberRunData &args, MVNDist *dist, string param_key)
{
    string filename = args.GetStringDefault(param_key, "modeldefault");
    if (filename != "modeldefault")
    {
        LOG << "VbInferenceTechnique::Loading " << param_key << " MVNDist from " << filename << endl;
        dist->LoadFromMatrix(filename);
        // Check the file we've loaded has the right number of parameters
        if (dist->GetSize() != m_num_params)
        {
            throw InvalidOptionValue(param_key, stringify(dist->GetSize()) + " parameters",
                "Model expected " + stringify(m_num_params) + " parameters");
        }
    }
}

void VariationalBayesInferenceTechnique::InitializeNoiseFromParam(FabberRunData &args, NoiseParams *dist,
    string param_key)
{
    string filename = args.GetStringDefault(param_key, "modeldefault");
    if (filename != "modeldefault")
    {
        LOG << "VbInferenceTechnique::Loading " << param_key << " distribution from " << filename << endl;
        dist->InputFromMVN(MVNDist(filename, m_log));
        // FIXME should there be checking of size here as well?
    }
}

void VariationalBayesInferenceTechnique::MakeInitialDistributions(FabberRunData &args)
{
    // Create initial prior and posterior distributions for model parameters
    initialFwdPrior = auto_ptr<MVNDist>(new MVNDist(m_num_params, m_log));
    initialFwdPosterior = new MVNDist(m_num_params, m_log);
    m_model->HardcodedInitialDists(*initialFwdPrior, *initialFwdPosterior);

    initialNoisePrior = noise->NewParams();
    initialNoisePosterior = noise->NewParams();
    noise->HardcodedInitialDists(*initialNoisePrior, *initialNoisePosterior);

    // The parameters fwd-initial-prior and fwd-initial-posterior can
    // be used to specify initial data for the prior and posterior
    // distributions.
    //
    // If unset, the default MVN distribution is used for each parameter
    // instead
    InitializeMVNFromParam(args, initialFwdPrior.get(), "fwd-initial-prior");
    InitializeMVNFromParam(args, initialFwdPosterior, "fwd-initial-posterior");

    // As above for the noise initial data
    InitializeNoiseFromParam(args, initialNoisePrior, "noise-initial-prior");
    InitializeNoiseFromParam(args, initialNoisePosterior, "noise-initial-posterior");
}

void VariationalBayesInferenceTechnique::GetPriorTypes(FabberRunData &args)
{
    vector<string> param_names;
    m_model->NameParams(param_names);

    m_prior_types.resize(0);
    for (int i = 0; i < m_num_params; i++)
    {
        PriorType prior = PriorType(i, param_names, args);
        LOG << prior;
        m_prior_types.push_back(prior);
        // FIXME HACK if a prior specifies a global precision, use it to override
        // the model default. Otherwise this is ignored in spatial VB. Note that this
        // only works because this method is called after MakeInitialDistributions. This
        // indicates a severe need to understand how priors are used in VB and SVB and
        // rationalize the whole mess.
        if (prior.m_prec > 0)
        {
            SymmetricMatrix prec = initialFwdPrior->GetPrecisions();
            prec(prior.m_idx + 1, prior.m_idx + 1) = prior.m_prec;
            initialFwdPrior->SetPrecisions(prec);
        }
    }
}

void VariationalBayesInferenceTechnique::Initialize(FwdModel *fwd_model, FabberRunData &args)
{
    // Call ancestor, which does most of the real work
    InferenceTechnique::Initialize(fwd_model, args);

    // Get noise model.
    noise = std::auto_ptr<NoiseModel>(NoiseModel::NewFromName(args.GetString("noise")));
    noise->Initialize(args);
    m_noise_params = noise->NumParams();
    LOG << "VariationalBayesInferenceTechnique::Noise has " << m_noise_params << " parameters" << endl;

    // Create initial prior and posterior distributions
    MakeInitialDistributions(args);

    // If we are resuming from a previous run, there will be a file containing a per-voxel
    // distribution of the model parameters, and possibly the noise as well. So we may
    // not need the initial posterior distributions we have created. We choose not to delete
    // them here as the memory involved is not large (they are not per voxel).
    m_continueFromFile = args.GetStringDefault("continue-from-mvn", "");
    paramFilename = args.GetStringDefault("continue-from-params", ""); // optional list of parameters in MVN

    // Get the spatial prior options for each parameter, if specified
    GetPriorTypes(args);

    // Create convergence-testing method:
    m_conv = ConvergenceDetector::NewFromName(args.GetStringDefault("convergence", "maxits"));
    m_conv->Initialize(args);

    // Figure out if F needs to be calculated every iteration
    m_printF = args.GetBool("print-free-energy");
    m_needF = m_conv->UseF() || m_printF;

    // Motion correction related setup - by default no motion correction
    Nmcstep = convertTo<int>(args.GetStringDefault("mcsteps", "0"));
}

void VariationalBayesInferenceTechnique::PassModelData(int voxel)
{
    // Pass in data, coords and supplemental data for this voxel
    ColumnVector y = m_origdata->Column(voxel);
    ColumnVector vcoords = m_coords->Column(voxel);
    if (m_suppdata->Ncols() > 0)
    {
        ColumnVector suppy = m_suppdata->Column(voxel);
        m_model->pass_in_data(y, suppy);
    }
    else
    {
        m_model->pass_in_data(y);
    }
    m_model->pass_in_coords(vcoords);
}

// continueFromDists = vector of MVNDist, one from each voxel, defining  parameter values from prev run
// modelpred = voxel data as predicted by the model at the end of each VB calculation. Used for motion correction
// data = copy of main voxel data
// resultMVNs = vector of MVNDist, one from each voxel, defining output parameter values after VB calculation
// resultFs = Free energy voxel data
// initialFwdPrior = Initial parameter distribution used for all voxels. Might have been loaded from file
//                   or obtained from model in Initialize()
// initialNoisePrior = Initial noise parameter distribution used for all voxels. Might have been loaded from file
//                     or obtained from model in Initialize()

void VariationalBayesInferenceTechnique::DoCalculations(FabberRunData &rundata)
{
    // Extract data (and the coords) from rundata for the (first) VB run
    // Number of rows is the size of the (time) series
    // Number of columns is the number of voxels
    m_origdata = &rundata.GetMainVoxelData();
    m_coords = &rundata.GetVoxelCoords();
    m_suppdata = &rundata.GetVoxelSuppData();
    m_nvoxels = m_origdata->Ncols();

#ifdef __FABBER_MOTION
    MCobj mcobj(rundata, 6); // hard coded DOF (future TODO item)
    // Copy the data for use in motion correction (FIXME why?)
    Matrix data(m_origdata->Nrows(), m_nvoxels);
    data = *m_origdata;
#endif //__FABBER_MOTION
    // Only call DoCalculations once
    assert(resultMVNs.empty());
    assert(resultFs.empty());

    // Initialize output data structures
    resultMVNs.resize(m_nvoxels, NULL);
    resultFs.resize(m_nvoxels, 9999); // 9999 is a garbage default value

    // Use this to store the model predictions in to pass to motion correction routine
    Matrix modelpred(m_origdata->Nrows(), m_nvoxels);

    //Indicates that we should continue from a previous run (e.g. after a motion correction step)
    bool continueFromPrevious = false;

    // If we're continuing from previous saved results, load the current
    // values of the parameters here: FIXME paramFilename is broken and ignored
    bool outputOnly = rundata.GetBool("output-only");
    if (m_continueFromFile != "")
    {
        LOG << "VbInferenceTechnique::Continuing from file " << m_continueFromFile << endl;
        InitMVNFromFile(m_continueFromFile, rundata, paramFilename);
        continueFromPrevious = true;
        if (outputOnly)
        {
            // Do no calculations - now we have set resultMVNs we can finish
            return;
        }
    }
    else if (outputOnly)
    {
        throw InvalidOptionValue("output-only", "", "Must also specify continue-from-mvn");
    }

    // Main loop over motion correction iterations and VB calculations
    for (int step = 0; step <= Nmcstep; step++)
    {
        if (step > 0)
            LOG << "VbInferenceTechnique::Motion correction step " << step << " of " << Nmcstep << endl;

        // Loop over voxels
        for (int voxel = 1; voxel <= m_nvoxels; voxel++)
        {
            PassModelData(voxel);

            // Apart from the data and the coords the main things we need here
            // are
            //   noiseVoxPrior		NoiseParams (prior)
            //   noiseVox			NoiseParams (posterior)
            //   fwdPrior			MVNDist
            //   fwdPosterior		MVNDist
            //
            // The posterior either comes from previous iteration (per-voxel)
            //   or from the file loaded from a previous run (per-voxel)
            //   or from initialNoisePosterior (same for all voxels)
            //     (which may have come from a file or the model!)
            // The prior comes from initialNoisePrior (same for all voxels)
            //     (which may have come from a file or the model!)
            //   or from an image prior (per-voxel) if specified

            const NoiseParams *noisePrior = initialNoisePrior;
            NoiseParams *noisePosterior = NULL;
            MVNDist fwdPrior(*initialFwdPrior);
            MVNDist fwdPosterior;

            for (int k = 0; k < m_num_params; k++)
            {
                m_prior_types[k].SetPrior(&fwdPrior, voxel);
            }

            if (continueFromPrevious)
            {
                // If we are continuing from a previous run still in memory, resultMVN will be populated
                // and contain model and noise parameters. Get the noise and model parameters from submatrices of
                // the result for this voxel
                noisePosterior = noise->NewParams();
                noisePosterior->InputFromMVN(
                    resultMVNs.at(voxel - 1)->GetSubmatrix(m_num_params + 1, m_num_params + m_noise_params));
                fwdPosterior = resultMVNs.at(voxel - 1)->GetSubmatrix(1, m_num_params);
            }
            else
            {
                // No previous run - just copy the standard initial posterior parameter values
                fwdPosterior = *initialFwdPosterior;
                noisePosterior = initialNoisePosterior->Clone();

                // The model may override the standard posterior on a per-voxel basis
                m_model->InitParams(fwdPosterior);
            }

            // Save our model parameters in case we need to revert later.
            // Note need to save prior in case ARD is being used
            NoiseParams *const noisePosteriorSave = noisePosterior->Clone();
            MVNDist fwdPosteriorSave(fwdPosterior);
            MVNDist fwdPriorSave(fwdPrior);

            // Give an indication of the progress through the voxels;
            rundata.Progress(voxel, m_nvoxels);
            double F = 1234.5678;

            // Create a linearized version of the model
            LinearizedFwdModel linear(m_model);

            // Setup for ARD (fwdmodel will decide if there is anything to be done)
            double Fard = 0;
            m_model->SetupARD(fwdPosterior, fwdPrior, Fard); // THIS USES ARD IN THE MODEL AND IS DEPRECEATED
            Fard = noise->SetupARD(m_model->ardindices, fwdPosterior, fwdPrior);

            try
            {
                linear.ReCentre(fwdPosterior.means);
                noise->Precalculate(*noisePosterior, *noisePrior, m_origdata->Column(voxel));
                m_conv->Reset();

                // START the VB updates and run through the relevant iterations (according to the convergence testing)
                int iteration = 0; //count the iterations
                do
                {
                    if (m_conv->NeedRevert()) //revert to previous solution if the convergence detector calls for it
                    {
                        *noisePosterior = *noisePosteriorSave;
                        fwdPosterior = fwdPosteriorSave;
                        fwdPrior = fwdPriorSave;
                        linear.ReCentre(fwdPosterior.means);
                    }

                    if (m_needF)
                    {
                        F = noise->CalcFreeEnergy(*noisePosterior, *noisePrior, fwdPosterior, fwdPrior, linear,
                            m_origdata->Column(voxel));
                        F = F + Fard;
                    }
                    if (m_printF)
                    {
                        LOG << "VbInferenceTechnique::Fbefore == " << F << endl;
                    }

                    // Save old values if called for
                    if (m_conv->NeedSave())
                    {
                        *noisePosteriorSave = *noisePosterior; // copy values, not pointer!
                        fwdPosteriorSave = fwdPosterior;
                        fwdPriorSave = fwdPrior;
                    }

                    // Do ARD updates (model will decide if there is anything to do here)
                    if (iteration > 0)
                    {
                        m_model->UpdateARD(fwdPosterior, fwdPrior, Fard); // THIS USES ARD IN THE MODEL AND IS DEPRECEATED
                        Fard = noise->UpdateARD(m_model->ardindices, fwdPosterior, fwdPrior);
                    }

                    // Theta update
                    noise->UpdateTheta(*noisePosterior, fwdPosterior, fwdPrior, linear, m_origdata->Column(voxel), NULL,
                        m_conv->LMalpha());

                    if (m_needF)
                    {
                        F = noise->CalcFreeEnergy(*noisePosterior, *noisePrior, fwdPosterior, fwdPrior, linear,
                            m_origdata->Column(voxel));
                        F = F + Fard;
                    }
                    if (m_printF)
                    {
                        LOG << "VbInferenceTechnique::Ftheta == " << F << endl;
                    }

                    // Alpha & Phi updates
                    noise->UpdateNoise(*noisePosterior, *noisePrior, fwdPosterior, linear, m_origdata->Column(voxel));

                    if (m_needF)
                    {
                        F = noise->CalcFreeEnergy(*noisePosterior, *noisePrior, fwdPosterior, fwdPrior, linear,
                            m_origdata->Column(voxel));
                        F = F + Fard;
                    }
                    if (m_printF)
                    {
                        LOG << "VbInferenceTechnique::Fphi == " << F << endl;
                    }

                    // Linearization update
                    // Update the linear model before doing Free energy calculation (and ready for next round of theta and phi updates)
                    linear.ReCentre(fwdPosterior.means);

                    if (m_needF)
                    {
                        F = noise->CalcFreeEnergy(*noisePosterior, *noisePrior, fwdPosterior, fwdPrior, linear,
                            m_origdata->Column(voxel));
                        F = F + Fard;
                    }
                    if (m_printF)
                    {
                        LOG << "VbInferenceTechnique::Fnoise == " << F << endl;
                    }

                    iteration++;
                } while (!m_conv->Test(F));

                // Revert to old values at last stage if required
                if (m_conv->NeedRevert())
                {
                    *noisePosterior = *noisePosteriorSave; // copy values, not pointers!
                    fwdPosterior = fwdPosteriorSave;
                    fwdPrior = fwdPriorSave;
                    linear.ReCentre(fwdPosterior.means); //just in case we go on to use this in motion correction
                }
            }
            catch (const overflow_error &e)
            {
                LOG_ERR("VbInferenceTechnique::Went infinite!  Reason:" << endl
                                                                        << "      " << e.what() << endl);
                //todo: write garbage or best guess to memory/file
                if (m_halt_bad_voxel)
                    throw;
                LOG_ERR("VbInferenceTechnique::Going on to the next voxel." << endl);
            }
            catch (Exception)
            {
                LOG_ERR("    NEWMAT Exception in this voxel:\n"
                    << Exception::what() << endl);
                if (m_halt_bad_voxel)
                    throw;
                LOG_ERR("VbInferenceTechnique::Going on to the next voxel." << endl);
            }
            catch (...)
            {
                LOG_ERR("VbInferenceTechnique::Other exception caught in main calculation loop!!\n");
                //<< "    Use --halt-on-bad-voxel for more details." << endl;
                if (m_halt_bad_voxel)
                    throw;
                LOG_ERR("VbInferenceTechnique::Going on to the next voxel" << endl);
            }

            // now write the results to resultMVNs
            try
            {
#if 0
				LOG << "VbInferenceTechnique::Final parameter estimates (" << fwdPosterior.means.Nrows() << "x"
				<< fwdPosterior.means.Ncols() << ") are: " << fwdPosterior.means.t() << endl;
				linear.DumpParameters(fwdPosterior.means, "      ");
#endif

                resultMVNs.at(voxel - 1) = new MVNDist(fwdPosterior, noisePosterior->OutputAsMVN());
                if (m_needF)
                {
                    resultFs.at(voxel - 1) = F;
                }
                // get the model prediction which is stored within the linearized forward model
                modelpred.Column(voxel) = linear.Offset();
            }
            catch (...)
            {
                // Even that can fail, due to results being singular
                LOG
                    << "VbInferenceTechnique::Can't give any sensible answer for this voxel; outputting zero +- identity\n";
                MVNDist *tmp = new MVNDist(m_log);
                tmp->SetSize(fwdPosterior.means.Nrows() + noisePosterior->OutputAsMVN().means.Nrows());
                tmp->SetCovariance(IdentityMatrix(tmp->means.Nrows()));
                resultMVNs.at(voxel - 1) = tmp;

                if (m_needF)
                    resultFs.at(voxel - 1) = F;
                modelpred.Column(voxel) = linear.Offset(); // get the model prediction which is stored within the linearized forward model
            }

            delete noisePosterior;
            noisePosterior = NULL;
            delete noisePosteriorSave;

        } //END of voxelwise updates

        //MOTION CORRECTION
        if (step < Nmcstep)
        { //dont do motion correction on the last run though as that would be a waste
#ifdef __FABBER_MOTION
            mcobj.run_mc(modelpred, data);
#endif //__FABBER_MOTION
        }

        continueFromPrevious = true; // We now take resultMVNs and use these as the starting point if we are to run again
    }                                // END of Steps that include motion correction and VB updates
}

void VariationalBayesInferenceTechnique::SaveResults(FabberRunData &rundata) const
{
    int nVoxels = resultMVNs.size();

    if (rundata.GetBool("save-noise-mean") | rundata.GetBool("save-noise-std"))
    {
        if (m_noise_params > 0)
        {
            LOG << "VariationalBayesInferenceTechnique::Writing noise" << endl;
            Matrix noiseMean, noiseStd;
            noiseMean.ReSize(m_noise_params, nVoxels);
            noiseStd.ReSize(m_noise_params, nVoxels);
            for (int vox = 1; vox <= nVoxels; vox++)
            {
                for (int i = 1; i <= m_noise_params; i++)
                {
                    noiseStd(i, vox) = sqrt(resultMVNs[vox - 1]->GetCovariance()(i + m_num_params, i + m_num_params));
                    noiseMean(i, vox) = resultMVNs[vox - 1]->means(i + m_num_params);
                }
            }
            // FIXME was this being saved before? Should it be?
            if (rundata.GetBool("save-noise-mean"))
                rundata.SaveVoxelData("noise_means", noiseMean);
            if (rundata.GetBool("save-noise-std"))
                rundata.SaveVoxelData("noise_stdevs", noiseStd);
        }
    }

    if (rundata.GetBool("save-mvn"))
    {
        if (resultMVNsWithoutPrior.size() > 0)
        {
            MVNDist::Save(resultMVNsWithoutPrior, "finalMVNwithoutPrior", rundata);
        }
    }

    // Save the Free Energy estimates
    if (rundata.GetBool("save-free-energy") && !resultFs.empty())
    {
        LOG << "VariationalBayesInferenceTechnique::Writing free energy" << endl;
        assert((int)resultFs.size() == nVoxels);
        Matrix freeEnergy;
        freeEnergy.ReSize(1, nVoxels);
        for (int vox = 1; vox <= nVoxels; vox++)
        {
            freeEnergy(1, vox) = resultFs.at(vox - 1);
        }
        rundata.SaveVoxelData("freeEnergy", freeEnergy);
    }
    else
    {
        LOG << "VariationalBayesInferenceTechnique::Free energy wasn't recorded, so no freeEnergy.nii.gz created.\n";
    }

    InferenceTechnique::SaveResults(rundata);
}

VariationalBayesInferenceTechnique::~VariationalBayesInferenceTechnique()
{
    delete m_conv;
    //delete initialFwdPrior;
    delete initialFwdPosterior;
    while (!resultMVNsWithoutPrior.empty())
    {
        delete resultMVNsWithoutPrior.back();
        resultMVNsWithoutPrior.pop_back();
    }
}
