/*  inference_vb.cc - VB inference technique class declarations

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_vb.h"
#include "convergence.h"

InferenceTechnique* VariationalBayesInferenceTechnique::NewInstance()
{
	return new VariationalBayesInferenceTechnique();
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
			throw Invalid_option("Size mismatch: model wants " + stringify(m_num_params) + ", loaded MVNDist has "
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
	MVNDist* loadPrior = new MVNDist(m_num_params);
	MVNDist* loadPosterior = new MVNDist(m_num_params);
	model->HardcodedInitialDists(*loadPrior, *loadPosterior);

	NoiseParams* loadNoisePrior = noise->NewParams();
	NoiseParams* loadNoisePosterior = noise->NewParams();
	noise->HardcodedInitialDists(*loadNoisePrior, *loadNoisePosterior);

	// The parameters fwd-initial-prior and fwd-initial-posterior can
	// be used to specify initial data for the prior and posterior
	// distributions.
	//
	// If unset, the default MVN distribution is used for each parameter
	// instead
	InitializeMVNFromParam(args, loadPrior, "fwd-initial-prior");
	InitializeMVNFromParam(args, loadPosterior, "fwd-initial-posterior");

	// As above for the noise initial data
	InitializeNoiseFromParam(args, loadNoisePrior, "noise-initial-prior");
	InitializeNoiseFromParam(args, loadNoisePosterior, "noise-initial-posterior");

	// Make these distributions constant: FIXME worth it?
	assert(initialFwdPrior == NULL);
	assert(initialFwdPosterior == NULL);
	assert(initialNoisePrior == NULL);
	assert(initialNoisePosterior == NULL);
	initialFwdPrior = loadPrior;
	initialFwdPosterior = loadPosterior;
	initialNoisePrior = loadNoisePrior;
	initialNoisePosterior = loadNoisePosterior;
	loadPrior = loadPosterior = NULL;
	loadNoisePrior = loadNoisePosterior = NULL; // now, only accessible as consts.
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
	vector<string> modnames;
	model->NameParams(modnames);
	imagepriorstr.resize(m_num_params);
	PriorsTypes.resize(m_num_params, ' ');
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
				char psp_type = convertTo<char> (args.GetString("PSP_byname" + stringify(current_psp) + "_type"));
				PriorsTypes[p] = psp_type;
				//record the index at which a PSP has been defined for use in spatialvb setup
				PSPidx.push_back(p);
				LOG << "VbInferenceTechnique::PSP_byname parameter " << param_name << " at entry " << p << ", type: "
						<< psp_type << ":" << PSPidx.size() << endl;

				// now read in file name for an image prior (if appropriate)
				if (psp_type == 'I')
				{
					imagepriorstr[p] = "PSP_byname" + stringify(current_psp) + "_image";
				}
			}
		}
		if (!found)
		{
			throw Invalid_option("ERROR: Prior specification by name, parameter " + param_name
					+ " does not exist in the model\n");
		}
	}
}

void VariationalBayesInferenceTechnique::Initialize(FwdModel* fwd_model, FabberRunData& args)
{
	Tracer_Plus tr("VariationalBayesInferenceTechnique::Initialize");

	// Call ancestor, which does most of the real work
	InferenceTechnique::Initialize(fwd_model, args);

	// Create initial prior and posterior distributions
	MakeInitialDistributions(args);

	// If we are resuming from a previous run, there will be a file containing a per-voxel
	// distribution of the model parameters, and possibly the noise as well. So we may
	// not need the initial posterior distributions we have created
	continueFromFile = args.GetStringDefault("continue-from-mvn", "");
	paramFilename = args.GetStringDefault("continue-from-params", ""); // optional list of parameters in MVN
	if (continueFromFile != "")
	{
		// Won't need these any more.  They don't hurt, but why leave them around?
		// They can only cause trouble (if used by mistake).
		delete initialFwdPosterior;
		initialFwdPosterior = NULL;

		if (!args.GetBool("continue-fwd-only"))
		{
			delete initialNoisePosterior;
			initialNoisePosterior = NULL;
		}
	}

	GetPriorTypes(args);

	// Create convergence-testing method:
	conv = ConvergenceDetector::NewFromName(args.GetStringDefault("convergence", "maxits"));
	conv->Initialize(args);

	// Fix the linearization centres? FIXME don't know what this is
	lockedLinearFile = args.GetStringDefault("locked-linear-from-mvn", "");

	// Figure out if F needs to be calculated every iteration
	printF = args.GetBool("print-free-energy");
	needF = conv->UseF() || printF;

	haltOnBadVoxel = !args.GetBool("allow-bad-voxels");
	if (haltOnBadVoxel)
	{
		LOG << "VbInferenceTechnique::Note: numerical errors in voxels will cause the program to halt.\n"
				<< "VbInferenceTechnique::Use --allow-bad-voxels (with caution!) to keep on calculating.\n";
	}
	else
	{
		LOG << "VbInferenceTechnique::Using --allow-bad-voxels: numerical errors in a voxel will\n"
				<< "VbInferenceTechnique::simply stop the calculation of that voxel.\n"
				<< "VbInferenceTechnique::Check log for 'Going on to the next voxel' messages.\n"
				<< "VbInferenceTechnique::Note that you should get very few (if any) exceptions like this;"
				<< "VbInferenceTechnique::they are probably due to bugs or a numerically unstable model.";
	}
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
	if (m_origdata->Nrows() != model->NumOutputs())
	{
		throw Invalid_option("Data length (" + stringify(m_origdata->Nrows())
				+ ") does not match model's output length (" + stringify(model->NumOutputs()) + ")!");
	}

	LoadImagePriors(allData);

#ifdef __FABBER_MOTION
	MCobj mcobj(allData,6); // hard coded DOF (future TODO item)
#endif //__FABBER_MOTION
	// Copy the data (FIXME why?)
	Matrix data(m_origdata->Nrows(), Nvoxels);
	data = *m_origdata;

	// Use this to store the model predictions in to pass to motion correction routine
	Matrix modelpred(model->NumOutputs(), Nvoxels);

	// Only call DoCalculations once
	assert(resultMVNs.empty());
	assert(resultFs.empty());

	// Initialize output data structures
	resultMVNs.resize(Nvoxels, NULL);
	resultFs.resize(Nvoxels, 9999); // 9999 is a garbage default value

	// If we're continuing from previous saved results, load the current
	// values of the parameters here: FIXME paramFilename is broken and ignored
	// FIXME could we just load directly into resultMVNs?

	bool continuingFromFile = (continueFromFile != "");
	if (continuingFromFile)
	{
		LOG << "VbInferenceTechnique::Continuing from file " << continueFromFile << endl;
		InitMVNFromFile(continueFromFile, allData, paramFilename);
	}

	// Don't know about this
	if (lockedLinearFile != "")
	{
		throw Invalid_option(
				"The option --locked-linear-from-mvn doesn't work with --method=vb yet, but should be pretty easy to implement.\n");
	}

	// Main loop over motion correction iterations and VB calculations
	bool continuefromprevious = false; //indicates that we should continue from a previous run (i.e. after a motion correction step)
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
			// A lot of the code below revolves around initializing these
			// things which may come from defaults, previous runs, files or image priors
			// Considerable improvement in readability would come from moving this
			// data initialization elsewhere and refactoring, once we have some
			// tests in place.

			// Create noise parameters for this voxel
			// The posterior either comes from previous iteration (per-voxel)
			//   or from the file loaded from a previous run (per-voxel)
			//   or from initialNoisePosterior (same for all voxels)
			//     (which may have come from a file or the model!)
			// The prior comes from initialNoisePrior
			//     (which may have come from a file or the model!)
			NoiseParams* noiseVox = NULL;

			if (continuefromprevious || continuingFromFile)
			{
				// If we are continuing from a previous run still in memory, resultMVN will be populated
				// and contain model and noise parameters. Get the noise parameters from the submatrix of
				// the result for this voxel
				noiseVox = noise->NewParams();
				noiseVox->InputFromMVN(resultMVNs.at(voxel - 1)->GetSubmatrix(m_num_params + 1, m_num_params
						+ m_noise_params));
			}
			else
			{
				// No previous run - just copy the standard initial noise values
				noiseVox = initialNoisePosterior->Clone();
			}

			const NoiseParams* noiseVoxPrior = initialNoisePrior;
			NoiseParams* const noiseVoxSave = noiseVox->Clone();

			// Give an indication of the progress through the voxels;
			allData.Progress(voxel, Nvoxels);
			double F = 1234.5678;

			// Determine the model parameters prior and posterior, in the same
			// way as the noise parameters were above
			MVNDist fwdPrior(*initialFwdPrior);
			MVNDist fwdPosterior;
			if (continuefromprevious || continuingFromFile)
			{
				// Similar to the extraction of noise posteriors above -
				// here we are using result from previous run still in
				// memory (presumably motion correction)
				fwdPosterior = resultMVNs.at(voxel - 1)->GetSubmatrix(1, m_num_params);
			}
			else
			{
				// No previous results - just use default
				assert(initialFwdPosterior != NULL);
				fwdPosterior = *initialFwdPosterior;
				// any voxelwise initialisation
				model->InitParams(fwdPosterior);
			}

			// Save our model parameters in case we need to revert later
			MVNDist fwdPosteriorSave(fwdPosterior);
			MVNDist fwdPriorSave(fwdPrior);

			// Create a linearized version of the model
			LinearizedFwdModel linear(model);

			// Setup for ARD (fwdmodel will decide if there is anything to be done)
			double Fard = 0;
			model->SetupARD(fwdPosterior, fwdPrior, Fard); // THIS USES ARD IN THE MODEL AND IS DEPRECEATED
			Fard = noise->SetupARD(model->ardindices, fwdPosterior, fwdPrior);

			// Image priors. This is a way of setting the mean of the
			// prior for each parameter on a per-voxel basis
			for (int k = 1; k <= m_num_params; k++)
			{
				if (PriorsTypes[k - 1] == 'I')
				{
					fwdPrior.means(k) = ImagePrior[k - 1](voxel);
				}
			}

			try
			{
				linear.ReCentre(fwdPosterior.means);
				noise->Precalculate(*noiseVox, *noiseVoxPrior, m_origdata->Column(voxel));
				conv->Reset();

				// START the VB updates and run through the relevant iterations (according to the convergence testing)
				int iteration = 0; //count the iterations
				do
				{
					if (conv-> NeedRevert()) //revert to previous solution if the convergence detector calls for it
					{
						*noiseVox = *noiseVoxSave; // copy values, not pointers!
						fwdPosterior = fwdPosteriorSave;
						fwdPrior = fwdPriorSave; // need to revert prior too (in case ARD is in place)
						linear.ReCentre(fwdPosterior.means);
					}

					if (needF)
					{
						F = noise->CalcFreeEnergy(*noiseVox, *noiseVoxPrior, fwdPosterior, fwdPrior, linear,
								m_origdata->Column(voxel));
						F = F + Fard;
					}
					if (printF)
						LOG << "VbInferenceTechnique::Fbefore == " << F << endl;

					// Save old values if called for
					if (conv->NeedSave())
					{
						*noiseVoxSave = *noiseVox; // copy values, not pointers!
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
					noise->UpdateTheta(*noiseVox, fwdPosterior, fwdPrior, linear, m_origdata->Column(voxel), NULL,
							conv->LMalpha());

					if (needF)
					{
						F = noise->CalcFreeEnergy(*noiseVox, *noiseVoxPrior, fwdPosterior, fwdPrior, linear,
								m_origdata->Column(voxel));
						F = F + Fard;
					}
					if (printF)
						LOG << "VbInferenceTechnique::Ftheta == " << F << endl;

					// Alpha & Phi updates
					noise->UpdateNoise(*noiseVox, *noiseVoxPrior, fwdPosterior, linear, m_origdata->Column(voxel));

					if (needF)
					{
						F = noise->CalcFreeEnergy(*noiseVox, *noiseVoxPrior, fwdPosterior, fwdPrior, linear,
								m_origdata->Column(voxel));
						F = F + Fard;
					}
					if (printF)
						LOG << "VbInferenceTechnique::Fphi == " << F << endl;

					// Test of NoiseModel cloning:
					// NoiseModel* tmp = noise; noise = tmp->Clone(); delete tmp;

					// Linearization update
					// Update the linear model before doing Free energy calculation (and ready for next round of theta and phi updates)
					linear.ReCentre(fwdPosterior.means);

					if (needF)
					{
						F = noise->CalcFreeEnergy(*noiseVox, *noiseVoxPrior, fwdPosterior, fwdPrior, linear,
								m_origdata->Column(voxel));
						F = F + Fard;
					}
					if (printF)
						LOG << "VbInferenceTechnique::Fnoise == " << F << endl;

					iteration++;
				} while (!conv->Test(F));

				// Revert to old values at last stage if required
				if (conv-> NeedRevert())
				{
					*noiseVox = *noiseVoxSave; // copy values, not pointers!
					fwdPosterior = fwdPosteriorSave;
					fwdPrior = fwdPriorSave;
					linear.ReCentre(fwdPosterior.means); //just in case we go on to use this in motion correction
				}
			} catch (const overflow_error& e)
			{
				LOG_ERR("VbInferenceTechnique::Went infinite!  Reason:" << endl
						<< "      " << e.what() << endl);
				//todo: write garbage or best guess to memory/file
				if (haltOnBadVoxel)
					throw;
				LOG_ERR("VbInferenceTechnique::Going on to the next voxel." << endl);
			} catch (Exception)
			{
				LOG_ERR("    NEWMAT Exception in this voxel:\n"
						<< Exception::what() << endl);
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

				resultMVNs.at(voxel - 1) = new MVNDist(fwdPosterior, noiseVox->OutputAsMVN());
				if (needF)
					resultFs.at(voxel - 1) = F;
				modelpred.Column(voxel) = linear.Offset(); // get the model prediction which is stored within the linearized forward model

			} catch (...)
			{
				// Even that can fail, due to results being singular
				LOG
						<< "VbInferenceTechnique::Can't give any sensible answer for this voxel; outputting zero +- identity\n";
				MVNDist* tmp = new MVNDist();
				tmp->SetSize(fwdPosterior.means.Nrows() + noiseVox->OutputAsMVN().means.Nrows());
				tmp->SetCovariance(IdentityMatrix(tmp->means.Nrows()));
				resultMVNs.at(voxel - 1) = tmp;

				if (needF)
					resultFs.at(voxel - 1) = F;
				modelpred.Column(voxel) = linear.Offset(); // get the model prediction which is stored within the linearized forward model
			}

			delete noiseVox;
			noiseVox = NULL;
			delete noiseVoxSave;
		} //END of voxelwise updates

		//MOTION CORRECTION
		if (step < Nmcstep)
		{ //dont do motion correction on the last run though as that would be a waste
#ifdef __FABBER_MOTION
			mcobj.run_mc(modelpred,data);
#endif //__FABBER_MOTION
		}

		continuefromprevious = true; // We now take resultMVNs and use these as the starting point if we are to run again
	} // END of Steps that include motion correction and VB updates

#if 0
	while (continueFromDists.size() > 0)
	{
		delete continueFromDists.back();
		continueFromDists.pop_back();
	}
#endif
}

VariationalBayesInferenceTechnique::~VariationalBayesInferenceTechnique()
{
	delete conv;
	delete initialFwdPrior;
	delete initialFwdPosterior;
}

