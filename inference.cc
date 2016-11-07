/*  inference.cc - General inference technique base class

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-20015 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference.h"
#include "fabber_mc.h"
#include "easylog.h"

#include "utils/tracer_plus.h"
#include "newmat.h"

#include "math.h"

using namespace std;
using namespace NEWMAT;
using Utilities::Tracer_Plus;

std::vector<std::string> InferenceTechnique::GetKnown()
{
	InferenceTechniqueFactory* factory = InferenceTechniqueFactory::GetInstance();
	return factory->GetNames();
}

InferenceTechnique* InferenceTechnique::NewFromName(const string& name)
{
	InferenceTechniqueFactory* factory = InferenceTechniqueFactory::GetInstance();
	InferenceTechnique* inf = factory->Create(name);
	if (inf == NULL)
	{
		throw Invalid_option("Unrecognized inference method: " + name);
	}
	return inf;
}

void InferenceTechnique::UsageFromName(const string& name, std::ostream &stream)
{
	stream << "Usage information for method: " << name << endl << endl;

	std::auto_ptr<InferenceTechnique> method(NewFromName(name));
	stream << method->GetDescription() << endl << endl << "Options: " << endl << endl;
	vector<OptionSpec> options;
	method->GetOptions(options);
	if (options.size() > 0)
	{
		for (vector<OptionSpec>::iterator iter = options.begin(); iter != options.end(); iter++)
		{
			stream << *iter << endl;
		}
	}
}

void InferenceTechnique::Initialize(FwdModel* fwd_model, FabberRunData& args)
{
	Tracer_Plus tr("InferenceTechnique::Initialize");

	// Set forward model.
	model = fwd_model;
	m_num_params = model->NumParams();
	LOG << "InferenceTechnique::Model has " << m_num_params << " parameters" << endl;

	saveModelFit = args.GetBool("save-model-fit");
	saveResiduals = args.GetBool("save-residuals");

	haltOnBadVoxel = !args.GetBool("allow-bad-voxels");
	if (haltOnBadVoxel)
	{
		LOG << "InferenceTechnique::Note: numerical errors in voxels will cause the program to halt.\n"
				<< "InferenceTechnique::Use --allow-bad-voxels (with caution!) to keep on calculating.\n";
	}
	else
	{
		LOG << "InferenceTechnique::Using --allow-bad-voxels: numerical errors in a voxel will\n"
				<< "InferenceTechnique::simply stop the calculation of that voxel.\n"
				<< "InferenceTechnique::Check log for 'Going on to the next voxel' messages.\n"
				<< "InferenceTechnique::Note that you should get very few (if any) exceptions like this;"
				<< "InferenceTechnique::they are probably due to bugs or a numerically unstable model.";
	}

	// Motion correction related setup - by default no motion correction
	Nmcstep = convertTo<int>(args.GetStringDefault("mcsteps", "0"));
}

void InferenceTechnique::SaveResults(FabberRunData& data) const
{
	Tracer_Plus tr("InferenceTechnique::SaveResults");
	LOG << "InferenceTechnique::Preparing to save results..." << endl;

	// Save the resultMVNs as two NIFTI files
	// Note: I should probably use a single NIFTI file with
	// NIFTI_INTENT_NORMAL -- but I can't find the detailed
	// documentation!  (Ordering for a multivariate norm).

	int nVoxels = resultMVNs.size();

	if (data.GetBool("save-mvn"))
	{
		MVNDist::Save(resultMVNs, "finalMVN", data);
		if (resultMVNsWithoutPrior.size() > 0)
		{
			assert(resultMVNsWithoutPrior.size() == (unsigned )nVoxels);
			MVNDist::Save(resultMVNsWithoutPrior, "finalMVNwithoutPrior", data);
		}
	}

	vector<string> paramNames;
	model->NameParams(paramNames);

#if 0
	LOG << "InferenceTechnique::Same information using DumpParameters:" << endl;
	ColumnVector indices(m_num_params);
	for (int i = 1; i <= indices.Nrows(); i++)
	indices(i) = i;

	model->DumpParameters(indices, "      ");
#endif

	// Create individual files for each parameter's mean and Z-stat

	if (data.GetBool("save-mean") | data.GetBool("save-std") | data.GetBool("save-zstat"))
	{

		LOG << "InferenceTechnique::Writing means..." << endl;
		for (unsigned i = 1; i <= paramNames.size(); i++)
		{
			Matrix paramMean, paramZstat, paramStd;
			paramMean.ReSize(1, nVoxels);
			paramZstat.ReSize(1, nVoxels);
			paramStd.ReSize(1, nVoxels);

			for (int vox = 1; vox <= nVoxels; vox++)
			{
				paramMean(1, vox) = resultMVNs[vox - 1]->means(i);
				float std = sqrt(resultMVNs[vox - 1]->GetCovariance()(i, i));
				paramZstat(1, vox) = paramMean(1, vox) / std;
				paramStd(1, vox) = std;
			}

			if (data.GetBool("save-mean"))
				data.SaveVoxelData("mean_" + paramNames.at(i - 1), paramMean);
			if (data.GetBool("save-zstat"))
				data.SaveVoxelData("zstat_" + paramNames.at(i - 1), paramZstat);
			if (data.GetBool("save-std"))
				data.SaveVoxelData("std_" + paramNames.at(i - 1), paramStd);
		}
	}

	// That's it! We've written our outputs to the "means" and "stdevs" output matrices.
	// Also save the noise parameters, just 'cuz.

	// We know how many noise parameters we have because it's the difference between
	// the size of the output matrix and the number of parameters in the model.
	if (data.GetBool("save-noise-mean") | data.GetBool("save-noise-std"))
	{

		const int nParams = paramNames.size();
		const int nNoise = resultMVNs[0]->means.Nrows() - paramNames.size();
		if (nNoise > 0)
		{
			LOG << "InferenceTechnique::Writing noise" << endl;
			Matrix noiseMean, noiseStd;
			noiseMean.ReSize(nNoise, nVoxels);
			noiseStd.ReSize(nNoise, nVoxels);
			for (int vox = 1; vox <= nVoxels; vox++)
			{
				for (int i = 1; i <= nNoise; i++)
				{
					noiseStd(i, vox) = sqrt(resultMVNs[vox - 1]->GetCovariance()(i + nParams, i + nParams));
					noiseMean(i, vox) = resultMVNs[vox - 1]->means(i + nParams);
				}
			}
			// FIXME was this being saved before? Should it be?
			if (data.GetBool("save-noise-mean"))
				data.SaveVoxelData("noise_means", noiseMean);
			if (data.GetBool("save-noise-std"))
				data.SaveVoxelData("noise_stdevs", noiseStd);
		}
	}

	// Save the Free Energy estimates
	if (data.GetBool("save-free-energy") && !resultFs.empty())
	{
		LOG << "InferenceTechnique::Writing free energy" << endl;
		assert((int )resultFs.size() == nVoxels);
		Matrix freeEnergy;
		freeEnergy.ReSize(1, nVoxels);
		for (int vox = 1; vox <= nVoxels; vox++)
		{
			freeEnergy(1, vox) = resultFs.at(vox - 1);
		}
		data.SaveVoxelData("freeEnergy", freeEnergy);
	}
	else
	{
		LOG << "InferenceTechnique::Free energy wasn't recorded, so no freeEnergy.nii.gz created.\n";
	}

	if (saveModelFit || saveResiduals)
	{
		LOG << "InferenceTechnique::Writing model fit/residuals..." << endl;
		// Produce the model fit and residual volumeserieses

		Matrix modelFit, residuals, datamtx, coords;
		datamtx = data.GetMainVoxelData(); // it is just possible that the model needs the data in its calculations
		coords = data.GetVoxelCoords();
		modelFit.ReSize(datamtx.Nrows(), nVoxels);
		ColumnVector tmp;
		for (int vox = 1; vox <= nVoxels; vox++)
		{
			// pass in stuff that the model might need
			ColumnVector y = datamtx.Column(vox);
			ColumnVector vcoords = coords.Column(vox);
			model->pass_in_data(y);
			model->pass_in_coords(vcoords);

			// do the evaluation
			model->Evaluate(resultMVNs.at(vox - 1)->means.Rows(1, m_num_params), tmp);
			modelFit.Column(vox) = tmp;
		}

		// FIXME no NIFTI write yet because not per-voxel data

		if (saveResiduals)
		{
			residuals = datamtx - modelFit;
			data.SaveVoxelData("residuals", residuals);
		}
		if (saveModelFit)
		{
			data.SaveVoxelData("modelfit", modelFit);
		}
	}

#if 0
	{
		LOG << "InferenceTechnique::Writing model variances..." << endl;
		Matrix modelStd;
		Matrix datamtx = data.GetMainVoxelData();
		modelStd.ReSize(datamtx.Nrows(), nVoxels);
		for (int vox = 1; vox <= nVoxels; vox++)
		{
			LinearizedFwdModel lin(model);
			lin.ReCentre(resultMVNs.at(vox - 1)->means.Rows(1, m_num_params));
			Matrix var = resultMVNs.at(vox - 1)->GetCovariance().SymSubMatrix(1, m_num_params);

			Matrix mvar = lin.Jacobian() * var * lin.Jacobian().t();
			ColumnVector tmp(datamtx.Nrows());
			for (int i = 1; i <= datamtx.Nrows(); i++)
			{
				tmp(i) = sqrt(mvar(i, i));
			}

			modelStd.Column(vox) = tmp;
		}
		data.SaveVoxelData("modelstd", modelStd);

	}
#endif

	LOG << "InferenceTechnique::Done writing results." << endl;
}

void InferenceTechnique::InitMVNFromFile(string continueFromFile, FabberRunData& allData, string paramFilename = "")
{
// Loads in a MVN to set it as inital values for inference
// can cope with the special scenario in which extra parameters have been added to the inference
	Tracer_Plus tr("InferenceTechnique::InitMVNFromFile");

	LOG << "InferenceTechnique::Merging supplied MVN with model intialization." << endl;

	if (paramFilename == "")
	{
		MVNDist::Load(resultMVNs, "continue-from-mvn", allData);
	}
	else
	{
		// load in parameters FIXME this will not work right now
		LOG << "InferenceTechnique::Parameters named in file" << endl;

		// Get a list of parameter names which were associated with the previous
		// run so we can merge the MVNs with the full parameter set
		string currparam;
		ifstream paramFile((paramFilename).c_str());
		if (!paramFile.good())
		{
			throw Invalid_option("Check filename of the parameter name file. ");
		}
		vector<string> paramNames;
		LOG << "InferenceTechnique::Parameters from previous run: " << endl;
		while (paramFile.good())
		{
			getline(paramFile, currparam);
			paramNames.push_back(currparam);
			LOG << currparam << endl;
		}
		paramNames.pop_back(); //remove final empty line assocaited with eof

		// get the parameters in the model
		vector<string> ModelparamNames;
		model->NameParams(ModelparamNames);
		LOG << "InferenceTechnique::Parameters named in model" << endl;
		for (int p = 0; p < m_num_params; p++)
		{
			LOG << ModelparamNames[p] << endl;
		}

		//load in the MVN
		vector<MVNDist*> MVNfile;
		MVNDist::Load(MVNfile, continueFromFile, allData);

		// Get defaults from the model. The prior is not used, the posterior is used
		// if we don't have a posterior for a parameter in the file
		MVNDist tempprior(m_num_params);
		MVNDist tempposterior(m_num_params);
		model->HardcodedInitialDists(tempprior, tempposterior);

		// go through the parameters in the model and either:
		// 1.) load the MVN from MVNfile if it is included, or
		// 2.) use the default value from the model

		// first work out where parameters in file MVN go in the model
		LOG << "InferenceTechnique::Matching parameters from file with model:" << endl;
		// Vector of bools for each parameter in the model. True if found in the file we're loading
		vector<bool> usefile(m_num_params, false);
		// If true, location of parameter in file (starting at 0)
		vector<int> oldloc(m_num_params, 0);
		// Vector of bools, one for each parameter in the file. True to flag matched to a file parameter
		vector<bool> hasmatched(m_num_params, false);
		for (int p = 0; p < m_num_params; p++)
		{
			usefile[p] = false;
			for (unsigned q = 0; q < paramNames.size(); q++)
			{
				if (ModelparamNames[p] == paramNames[q])
				{
					usefile[p] = true;
					oldloc[p] = q;
					hasmatched[q] = true;
					LOG << ModelparamNames[p] << ": Matched with file" << endl;
				}
			}
			if (!usefile[p])
			{
				LOG << ModelparamNames[p] << ": Not matched, set from model default" << endl;
			}
		}

		//Make a note of any parameters in the file that were not matched to params in the model
		for (unsigned int q = 0; q < paramNames.size(); q++)
		{
			if (!hasmatched[q])
			{
				LOG_ERR(paramNames[q] + ": Not matched!");
			}
		}

		// for (int a=0; a<usefile.size(); a++) {
		//       cout << usefile[a] << "  " << oldloc[a] << endl;
		//     }

		ColumnVector modelmeans = tempposterior.means;
		ColumnVector newmeans = modelmeans;
		SymmetricMatrix modelcov = tempposterior.GetCovariance();
		SymmetricMatrix newcov = modelcov;

		// Number of forward model and noise params in the FILE we are loading
		int n_file_params = paramNames.size();
		int n_file_noiseparams = MVNfile[1]->means.Nrows() - n_file_params;
		int nvox = MVNfile.size();

		// For every voxel we need to add an MVNDist to resultMVNs.
		for (int v = 0; v < nvox; v++)
		{

			// Get the file data
			MVNDist fwddist = MVNfile[v]->GetSubmatrix(1, n_file_params);
			MVNDist noisedist = MVNfile[v]->GetSubmatrix(n_file_params + 1, n_file_params + n_file_noiseparams);

			for (unsigned int p = 0; p < ModelparamNames.size(); p++)
			{
				// deal with the means
				if (usefile[p])
				{
					newmeans(p + 1) = fwddist.means(oldloc[p] + 1);
				}
			}
			MVNDist newfwd(m_num_params);
			newfwd.means = newmeans;

			//deal with the covariances
			SymmetricMatrix filecov = fwddist.GetCovariance();
			for (unsigned int p = 0; p < ModelparamNames.size(); p++)
			{
				for (unsigned int q = 0; q <= p; q++)
				{
					if (usefile[p])
					{
						if (usefile[q])
						{
							newcov(p + 1, q + 1) = filecov(oldloc[p] + 1, oldloc[q] + 1);
						}
					}
				}
			}
			newfwd.SetCovariance(newcov);

			// Note that this assumes we have the same number of noise parameters in the
			// file as we want in the calculation, not sure at moment if this is always
			// true or not.
			MVNDist *distout = new MVNDist(newfwd, noisedist);
			resultMVNs.push_back(distout);
		}
	}
}

InferenceTechnique::~InferenceTechnique()
{
	while (!resultMVNs.empty())
	{
		delete resultMVNs.back();
		resultMVNs.pop_back();
	}
	while (!resultMVNsWithoutPrior.empty())
	{
		delete resultMVNsWithoutPrior.back();
		resultMVNsWithoutPrior.pop_back();
	}
}
