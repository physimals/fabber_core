/*  inference_vb.cc - VB inference technique class declarations

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_vb.h"

#include "convergence.h"

#include "utils/tracer_plus.h"
#include "newmatio.h"

using Utilities::Tracer_Plus;

static OptionSpec OPTIONS[] =
		{
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
				{ "continue-from-mvn", OPT_FILE, "Continue previous run from output MVN files", OPT_NONREQ, "" },
				{ "fwd-initial-prior", OPT_FILE, "MVN file containing initial model prior", OPT_NONREQ, "" },
				{ "fwd-initial-posterior", OPT_FILE, "MVN file containing initial model posterior", OPT_NONREQ, "" },
				{ "noise-initial-prior", OPT_FILE, "MVN file containing initial noise prior", OPT_NONREQ, "" },
				{ "noise-initial-posterior", OPT_FILE, "MVN file containing initial noise posterior", OPT_NONREQ, "" },
				{ "noise-pattern", OPT_STR,
						"repeating pattern of noise variances for each point (e.g. 12 gives odd and even data points different noise variances)",
						OPT_NONREQ, "1" },
				{ "PSP_byname<n>", OPT_STR, "Name of model parameter to use image prior", OPT_NONREQ, "" },
				{ "PSP_byname<n>_type", OPT_STR, "Type of image prior to use fo parameter <n> - I=image prior",
						OPT_NONREQ, "" },
				{ "PSP_byname<n>_image", OPT_FILE, "File containing image prior for parameter <n>", OPT_NONREQ, "" },
				{ "PSP_byname<n>_prec", OPT_FILE, "Precision to apply to image prior for parameter <n>", OPT_NONREQ, "" },
				{ "locked-linear-from-mvn", OPT_FILE, "MVN file containing fixed centres for linearization", OPT_NONREQ,
						"" },
				{ "allow-bad-voxels", OPT_BOOL, "Continue if numerical error found in a voxel, rather than stopping",
						OPT_NONREQ, "" },
				{ "ar1-cross-terms", OPT_STR, "For AR1 noise, type of cross-linking (dual, same or none)", OPT_NONREQ,
						"dual" },
				{ "" }, };

InferenceTechnique* VariationalBayesInferenceTechnique::NewInstance()
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

void VariationalBayesInferenceTechnique::InitializeMVNFromParam(FabberRunData& args, MVNDist *dist, string param_key)
{
	string filename = args.GetStringDefault(param_key, "modeldefault");
	if (filename != "modeldefault")
	{
		LOG << "VbInferenceTechnique::Loading " << param_key << " MVNDist from " << filename << endl;
		dist->Load(filename);
		// Check the file we've loaded has the right number of parameters
		if (dist->GetSize() != m_num_params)
		{
			throw Invalid_option(
					"Size mismatch: model wants " + stringify(m_num_params) + ", loaded MVNDist has "
							+ stringify(dist->GetSize()));
		}
	}
}

void VariationalBayesInferenceTechnique::InitializeNoiseFromParam(FabberRunData& args, NoiseParams *dist,
		string param_key)
{
	string filename = args.GetStringDefault(param_key, "modeldefault");
	if (filename != "modeldefault")
	{
		LOG << "VbInferenceTechnique::Loading " << param_key << " distribution from " << filename << endl;
		dist->InputFromMVN(MVNDist(filename));
		// FIXME should there be checking of size here as well?
	}
}

void VariationalBayesInferenceTechnique::MakeInitialDistributions(FabberRunData& args)
{
	// Create initial prior and posterior distributions for model parameters
	initialFwdPrior = auto_ptr < MVNDist > (new MVNDist(m_num_params));
	initialFwdPosterior = new MVNDist(m_num_params);
	model->HardcodedInitialDists(*initialFwdPrior, *initialFwdPosterior);

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

void VariationalBayesInferenceTechnique::GetPriorTypes(FabberRunData& args)
{
	// Get the prior types and sources for all parameters.
	//
	// Parameters used are:
	//   PSP_byname<n> = name of parameter
	//   PSP_byname<n>_type = 'I' for an image prior
	//   PSP_byname<n>_image = Filename/data key name for image prior data
	//
	// PriorsTypes is indexed by parameter number and contains the prior type.
	// imagepriorstr similarly contains the image file name/key

	// Get names of model parameters
	vector < string > modnames;
	model->NameParams(modnames);
	imagepriorstr.resize(m_num_params);
	PriorsTypes.resize(m_num_params, ' ');
	PriorsPrec.resize(m_num_params, 0);
	int current_psp = 0;
	while (true)
	{
		current_psp++;
		string param_name = args.GetStringDefault("PSP_byname" + stringify(current_psp), "stop!");
		if (param_name == "stop!")
			break; //no more spriors have been specified

		// Compare name to those in list of model names
		bool found = false;
		for (int p = 0; p < m_num_params; p++)
		{
			if (param_name == modnames[p])
			{
				found = true;
				char psp_type = convertTo<char>(args.GetString("PSP_byname" + stringify(current_psp) + "_type"));
				PriorsTypes[p] = psp_type;

				// A precision is optional. 0 means not set
				double prec = convertTo<double>(
						args.GetStringDefault("PSP_byname" + stringify(current_psp) + "_prec", "0"));
				PriorsPrec[p] = prec;

				// Record the index at which a PSP has been defined for use in spatialvb setup
				PSPidx.push_back(p);
				LOG << "VbInferenceTechnique::PSP_byname parameter " << param_name << " at entry " << p << ", type: "
						<< psp_type << ":" << PSPidx.size() << endl;

				// Record the data key for an image prior (if appropriate)
				if (psp_type == 'I')
				{
					imagepriorstr[p] = "PSP_byname" + stringify(current_psp) + "_image";
				}
			}
		}
		if (!found)
		{
			throw Invalid_option(
					"ERROR: Prior specification by name, parameter " + param_name + " does not exist in the model\n");
		}
	}
}

void VariationalBayesInferenceTechnique::Initialize(FwdModel* fwd_model, FabberRunData& args)
{
	Tracer_Plus tr("VariationalBayesInferenceTechnique::Initialize");

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

	// Fix the linearization centres? FIXME don't know what this is
	lockedLinearFile = args.GetStringDefault("locked-linear-from-mvn", "");
	if (lockedLinearFile != "")
	{
		throw Invalid_option(
				"The option --locked-linear-from-mvn doesn't work with --method=vb yet, but should be pretty easy to implement.\n");
	}

	// Figure out if F needs to be calculated every iteration
	m_printF = args.GetBool("print-free-energy");
	m_needF = m_conv->UseF() || m_printF;

}

void VariationalBayesInferenceTechnique::LoadImagePriors(FabberRunData &allData)
{
	// Load image priors. These are prior values for the parameters which
	// vary from voxel to voxel. Any selection of parameters can have image
	// priors.
	//
	// RowVector used as general convention is that voxel data is
	// stored as a one column per voxel.
	ImagePrior.resize(m_num_params);
	for (int k = 1; k <= m_num_params; k++)
	{
		LOG << "VbInferenceTechnique::Prior type for param " << k << " is " << PriorsTypes[k - 1] << endl;
		if (PriorsTypes[k - 1] == 'I')
		{
			string fname = imagepriorstr[k - 1];
			LOG << "VbInferenceTechnique::Reading Image prior (" << k << "): " << fname << endl;
			ImagePrior[k - 1] = allData.GetVoxelData(fname).AsRow();
		}
	}
}

void VariationalBayesInferenceTechnique::PassModelData(int voxel)
{
	// Pass in data, coords and supplemental data for this voxel
	ColumnVector y = m_origdata->Column(voxel);
	ColumnVector vcoords = m_coords->Column(voxel);
	if (m_suppdata->Ncols() > 0)
	{
		ColumnVector suppy = m_suppdata->Column(voxel);
		model->pass_in_data(y, suppy);
	}
	else
	{
		model->pass_in_data(y);
	}
	model->pass_in_coords(vcoords);
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

void VariationalBayesInferenceTechnique::DoCalculations(FabberRunData& allData)
{
	Tracer_Plus tr("VariationalBayesInferenceTechnique::DoCalculations");

	// extract data (and the coords) from allData for the (first) VB run
	// Rows are volumes
	// Columns are (time) series
	// num Rows is size of (time) series
	// num Cols is size of volumes
	m_origdata = &allData.GetMainVoxelData();
	m_coords = &allData.GetVoxelCoords();
	m_suppdata = &allData.GetVoxelSuppData();
	int Nvoxels = m_origdata->Ncols();

	// pass in some (dummy) data/coords here just in case the model relies upon it
	// use the first voxel values as our dummies FIXME this shouldn't really be
	// necessary, need to find way for model to know about the data beforehand.
	PassModelData(1);

	LoadImagePriors(allData);

#ifdef __FABBER_MOTION
	MCobj mcobj(allData,6); // hard coded DOF (future TODO item)
#endif //__FABBER_MOTION
	// Copy the data (FIXME why?)
	Matrix data(m_origdata->Nrows(), Nvoxels);
	data = *m_origdata;

	// Use this to store the model predictions in to pass to motion correction routine
	Matrix modelpred(m_origdata->Nrows(), Nvoxels);

	// Only call DoCalculations once
	assert(resultMVNs.empty());
	assert(resultFs.empty());

	// Initialize output data structures
	resultMVNs.resize(Nvoxels, NULL);
	resultFs.resize(Nvoxels, 9999); // 9999 is a garbage default value

	//Indicates that we should continue from a previous run (i.e. after a motion correction step)
	bool continueFromPrevious = false;

	// If we're continuing from previous saved results, load the current
	// values of the parameters here: FIXME paramFilename is broken and ignored
	if (m_continueFromFile != "")
	{
		LOG << "VbInferenceTechnique::Continuing from file " << m_continueFromFile << endl;
		InitMVNFromFile(m_continueFromFile, allData, paramFilename);
		continueFromPrevious = true;
	}

	// Main loop over motion correction iterations and VB calculations
	for (int step = 0; step <= Nmcstep; step++)
	{
		if (step > 0)
			LOG << "VbInferenceTechnique::Motion correction step " << step << " of " << Nmcstep << endl;

		// Loop over voxels
		for (int voxel = 1; voxel <= Nvoxels; voxel++)
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

			const NoiseParams* noisePrior = initialNoisePrior;
			NoiseParams* noisePosterior = NULL;
			MVNDist fwdPrior(*initialFwdPrior);
			MVNDist fwdPosterior;

			// We may have an image priors. This is a way of setting the mean of the
			// prior for each parameter on a per-voxel basis. We can also optionally
			// set the precision on an image prior parameter
			for (int k = 1; k <= m_num_params; k++)
			{
				if (PriorsTypes[k - 1] == 'I')
				{
					fwdPrior.means(k) = ImagePrior[k - 1](voxel);
					SymmetricMatrix prec = fwdPrior.GetPrecisions();
					bool changed = false;
					if (PriorsPrec[k - 1] != 0)
					{
						prec(k, k) = PriorsPrec[k - 1];
						changed = true;
					}
					LOG << prec << "," << prec.Nrows() << "," << prec.Ncols() << endl;
					if (changed)
						fwdPrior.SetPrecisions(prec);
				}
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
				model->InitParams(fwdPosterior);
			}

			// Save our model parameters in case we need to revert later.
			// Note need to save prior in case ARD is being used
			NoiseParams* const noisePosteriorSave = noisePosterior->Clone();
			MVNDist fwdPosteriorSave(fwdPosterior);
			MVNDist fwdPriorSave(fwdPrior);

			// Give an indication of the progress through the voxels;
			allData.Progress(voxel, Nvoxels);
			double F = 1234.5678;

			// Create a linearized version of the model
			LinearizedFwdModel linear(model);

			// Setup for ARD (fwdmodel will decide if there is anything to be done)
			double Fard = 0;
			model->SetupARD(fwdPosterior, fwdPrior, Fard); // THIS USES ARD IN THE MODEL AND IS DEPRECEATED
			Fard = noise->SetupARD(model->ardindices, fwdPosterior, fwdPrior);

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
						model->UpdateARD(fwdPosterior, fwdPrior, Fard); // THIS USES ARD IN THE MODEL AND IS DEPRECEATED
						Fard = noise->UpdateARD(model->ardindices, fwdPosterior, fwdPrior);
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
			} catch (const overflow_error& e)
			{
				LOG_ERR("VbInferenceTechnique::Went infinite!  Reason:" << endl << "      " << e.what() << endl);
				//todo: write garbage or best guess to memory/file
				if (haltOnBadVoxel)
					throw;
				LOG_ERR("VbInferenceTechnique::Going on to the next voxel." << endl);
			} catch (Exception)
			{
				LOG_ERR("    NEWMAT Exception in this voxel:\n" << Exception::what() << endl);
				if (haltOnBadVoxel)
					throw;
				LOG_ERR("VbInferenceTechnique::Going on to the next voxel." << endl);
			} catch (...)
			{
				LOG_ERR("VbInferenceTechnique::Other exception caught in main calculation loop!!\n");
				//<< "    Use --halt-on-bad-voxel for more details." << endl;
				if (haltOnBadVoxel)
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

			} catch (...)
			{
				// Even that can fail, due to results being singular
				LOG
						<< "VbInferenceTechnique::Can't give any sensible answer for this voxel; outputting zero +- identity\n";
				MVNDist* tmp = new MVNDist();
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
		mcobj.run_mc(modelpred,data);
#endif //__FABBER_MOTION
		}

		continueFromPrevious = true; // We now take resultMVNs and use these as the starting point if we are to run again
	} // END of Steps that include motion correction and VB updates
}

VariationalBayesInferenceTechnique::~VariationalBayesInferenceTechnique()
{
	delete m_conv;
	//delete initialFwdPrior;
	delete initialFwdPosterior;
}

