/*  inference.cc - General inference technique base class

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference.h"
#include "newimage/newimageall.h"
 
using namespace NEWIMAGE;
using namespace std;
using namespace MISCMATHS;

void InferenceTechnique::Setup(ArgsType& args)
{
  Tracer_Plus tr("InferenceTechnique::Setup");

  // Pick models
  model = FwdModel::NewFromName(args.Read("model"), args);
  assert( model->NumParams() > 0 );
  LOG_ERR("    Forward Model version:\n      " 
	  << model->ModelVersion() << endl);

  noise = NoiseModel::NewFromName(args.Read("noise"), args);
//  noise->LoadPrior(args.ReadWithDefault("noise-prior","hardcoded"));
//  noise->Dump("  ");

  saveModelFit = args.ReadBool("save-model-fit");
  saveResiduals = args.ReadBool("save-residuals");

  // Motion correction related setup
  Nmcstep = convertTo<int>(args.ReadWithDefault("mcsteps","0")); //by default no motion correction
}



void InferenceTechnique::SaveResults(const DataSet& data) const
{
  Tracer_Plus tr("InferenceTechnique::SaveResults");
    LOG << "    Preparing to save results..." << endl;

  
#ifdef __FABBER_LIBRARYONLY
    throw Logic_error("SaveResults shouldn't be called in fabber_library mode");
#else // __FABBER_LIBRARYONLY

    // Save the resultMVNs as two NIFTI files
    // Note: I should probably use a single NIFTI file with
    // NIFTI_INTENT_NORMAL -- but I can't find the detailed 
    // documentation!  (Ordering for a multivariate norm).

    const volume<float>& mask  = data.GetMask();
    int nVoxels = resultMVNs.size();

    cout << "Saving!\n";
    MVNDist::Save(resultMVNs, outputDir + "/finalMVN", mask);

    if (resultMVNsWithoutPrior.size() > 0)
      {
	assert(resultMVNsWithoutPrior.size() == (unsigned)nVoxels);
	MVNDist::Save(resultMVNsWithoutPrior, outputDir + "/finalMVNwithoutPrior", mask);
      }

    /* Some validation code -- checked, Save then Load 
       produced identical results (to single precision)
       cout << "Creating!\n";    
       vector<MVNDist*> test(resultMVNs.size(), NULL);
       cout << "Loading!\n";    
       MVNDist::Load(test, outputDir + "/finalMVN", mask);

       assert(test[0] != NULL);

       cout << "Verifying MVNDists are identical!!!";    
       // won't be identical because they're written as floats.
       for (unsigned i = 1; i <= test.size(); i++)
       {
       cout << i << endl;
       test.at(i-1); resultMVNs.at(i-1);
       cout << 'a'<< endl;
       assert(resultMVNs[i-1] != NULL);
       assert(test[i-1] != NULL);
       cout << resultMVNs[i-1]->means.t() << test[i-1]->means.t();
       //        assert(resultMVNs[i-1]->means == test[i-1]->means);
       cout << 'b' << endl;
       cout << resultMVNs[i-1]->GetCovariance();
       cout << test[i-1]->GetCovariance();
       //        assert(resultMVNs[i-1]->GetCovariance() == test[i-1]->GetCovariance());
       cout << 'c' << endl;
       }
    */

    // Write the parameter names into paramnames.txt
    
    LOG << "    Writing paramnames.txt..." << endl;
    ofstream paramFile((outputDir + "/paramnames.txt").c_str());
    vector<string> paramNames;
    model->NameParams(paramNames);
    for (unsigned i = 0; i < paramNames.size(); i++)
    {
        LOG << "      " << paramNames[i] << endl;
        paramFile << paramNames[i] << endl;
    }
    paramFile.close();
    
    LOG << "    Same information using DumpParameters:" << endl;
    ColumnVector indices(model->NumParams());
    for (int i = 1; i <= indices.Nrows(); i++)
        indices(i) = i;
    model->DumpParameters(indices, "      ");

    // Create individual files for each parameter's mean and Z-stat

    if (!EasyOptions::UsingMatrixIO())
    {
      for (unsigned i = 1; i <= paramNames.size(); i++)
      {
        Matrix paramMean, paramZstat;
    	paramMean.ReSize(1, nVoxels);
	paramZstat.ReSize(1, nVoxels);

        for (int vox = 1; vox <= nVoxels; vox++)
        {
	    paramMean(1,vox) = resultMVNs[vox-1]->means(i);
            paramZstat(1,vox) =
              paramMean(1,vox) / 
              sqrt(resultMVNs[vox-1]->GetCovariance()(i,i));
        }
    	LOG << "    Writing means..." << endl;

	volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),1);
	output.setmatrix(paramMean,mask);
	output.set_intent(NIFTI_INTENT_NONE,0,0,0);
	output.setDisplayMaximumMinimum(output.max(),output.min());
	save_volume4D(output,outputDir + "/mean_" + paramNames.at(i-1));

	output.setmatrix(paramZstat,mask);
        output.set_intent(NIFTI_INTENT_ZSCORE,0,0,0);
	output.setDisplayMaximumMinimum(output.max(),output.min());
	save_volume4D(output,outputDir + "/zstat_" + paramNames.at(i-1));
      }
    }        
    else
    {
	Matrix& paramMean = EasyOptions::OutMatrix("<means>"); // Creates matrix
	Matrix& paramStd = EasyOptions::OutMatrix("<stdevs>");
	const int nParams = paramNames.size();
	paramMean.ReSize(nParams, nVoxels);
	paramStd.ReSize(nParams, nVoxels);
	for (int vox = 1; vox <= nVoxels; vox++)
	{
            for (int i = 1; i <= nParams; i++)
	    {
              	paramStd(i,vox) = sqrt(resultMVNs[vox-1]->GetCovariance()(i,i));
	    	paramMean(i,vox) = resultMVNs[vox-1]->means(i);
	    }
	}
	// That's it! We've written our outputs to the "means" and "stdevs" output matrices.
	// Also save the noise parameters, just 'cuz.

	Matrix& noiseMean = EasyOptions::OutMatrix("<noise_means>"); // Creates matrix
	Matrix& noiseStd = EasyOptions::OutMatrix("<noise_stdevs>");
	const int nNoise = resultMVNs[0]->means.Nrows() - paramNames.size();
	noiseMean.ReSize(nNoise, nVoxels);
	noiseStd.ReSize(nNoise, nVoxels);
	for (int vox = 1; vox <= nVoxels; vox++)
	{
            for (int i = 1; i <= nNoise; i++)
	    {
              	noiseStd(i,vox) = sqrt(resultMVNs[vox-1]->GetCovariance()(i+nParams,i+nParams));
	    	noiseMean(i,vox) = resultMVNs[vox-1]->means(i+nParams);
	    }
	}
    }

    // Save the Free Energy estimates
    if (!resultFs.empty())
      {
	assert((int)resultFs.size() == nVoxels);
	Matrix freeEnergy;
	freeEnergy.ReSize(1, nVoxels);
	for (int vox = 1; vox <= nVoxels; vox++)
	  {
	    freeEnergy(1,vox) = resultFs.at(vox-1);
	  }
	
	if (EasyOptions::UsingMatrixIO())
	{
	  EasyOptions::OutMatrix("<freeEnergy>") = freeEnergy;
        }
        else
        {
	volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),1);
	output.setmatrix(freeEnergy,mask);
	output.set_intent(NIFTI_INTENT_NONE,0,0,0);
	output.setDisplayMaximumMinimum(output.max(),output.min());
	save_volume4D(output,outputDir + "/freeEnergy");
        }
      }
    else
      {
	LOG_ERR("Free energy wasn't recorded, so no freeEnergy.nii.gz created.\n");
      }
    
    if (saveModelFit || saveResiduals)
      {
        LOG << "    Writing model fit/residuals..." << endl;
        // Produce the model fit and residual volumeserieses
	
        Matrix modelFit, residuals, datamtx, coords;
        modelFit.ReSize(model->NumOutputs(), nVoxels);
	datamtx = data.GetVoxelData(); // it is just possible that the model needs the data in its calculations
	coords = data.GetVoxelCoords();
	ColumnVector tmp;
        for (int vox = 1; vox <= nVoxels; vox++)
        {
	  // pass in stuff that the model might need
	  ColumnVector y = datamtx.Column(vox);
	  ColumnVector vcoords = coords.Column(vox);
	  model->pass_in_data( y );
	  model->pass_in_coords(vcoords);

	  // do the evaluation
	  model->Evaluate(resultMVNs.at(vox-1)->means.Rows(1,model->NumParams()), tmp);
	  modelFit.Column(vox) = tmp;
        }

	volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),model->NumOutputs());
	
        if (saveResiduals)
        {
	  residuals = datamtx - modelFit;

         if (EasyOptions::UsingMatrixIO())
         {
          EasyOptions::OutMatrix("<residuals>") = residuals;
         }
	 else
	 { 
	  output.setmatrix(residuals,mask);
	  output.set_intent(NIFTI_INTENT_NONE,0,0,0);
	  output.setDisplayMaximumMinimum(output.max(),output.min());
	  save_volume4D(output,outputDir + "/residuals");
	 }
        }
        if (saveModelFit)
        {
	 if (EasyOptions::UsingMatrixIO())
	 {
	  EasyOptions::OutMatrix("<modelfit>") = modelFit;
	 }
	 else
	 {
 	  output.setmatrix(modelFit,mask);
	  output.set_intent(NIFTI_INTENT_NONE,0,0,0);
	  output.setDisplayMaximumMinimum(output.max(),output.min());
	  save_volume4D(output,outputDir + "/modelfit");
	 }
        }

    }

#endif // __FABBER_LIBRARYONLY
    LOG << "    Done writing results." << endl;
}

void InferenceTechnique::InitMVNFromFile(vector<MVNDist*>& continueFromDists,string continueFromFile, const DataSet& allData, string paramFilename="") {
#ifdef __FABBER_LIBRARYONLY
  throw Logic_error("Should not be called when compiled without NEWIMAGE support");
#else
  // Loads in a MVN to set it as inital values for inference
  // can cope with the special scenario in which extra parameters have been added to the inference
  Tracer_Plus tr("InferenceTechnique::InitMVNFromFile");

  LOG << "Merging supplied MVN with model intialization." << endl;

  if (paramFilename == "") {
    MVNDist::Load(continueFromDists, continueFromFile, allData.GetMask());
  }
  else {
    // load in parameters
    LOG << "Parameters named in file" << endl;
    string currparam;
    ifstream paramFile((paramFilename).c_str());
    if (!paramFile.good()) {
      throw Invalid_option("Check filename of the parameter name file. ");
    }
    vector<string> paramNames;
    while (paramFile.good())
      {
	getline(paramFile,currparam);
	paramNames.push_back(currparam);
	LOG << currparam << endl;
      }
    paramNames.pop_back(); //remove final empty line assocaited with eof
    
    // get the parameters in the model
    vector<string> ModelparamNames;
    model->NameParams(ModelparamNames);
    int nmodparams = model->NumParams();
    LOG << "Parameters named in model" << endl;
    for (int p=0; p<nmodparams; p++) {
      LOG << ModelparamNames[p] << endl;
    }

    //load in the MVN
    vector<MVNDist*> MVNfile;
    MVNDist::Load(MVNfile, continueFromFile, allData.GetMask());

    // Get deafults from the model
    
    MVNDist tempprior(nmodparams);
    MVNDist tempposterior(nmodparams);
    model->HardcodedInitialDists(tempprior,tempposterior);


    // go through the parameters in the model and either:
    // 1.) load the MVN from MVNfile if it is included, or
    // 2.) use the default value from the model

    // first work out where parameters in file MVN go in the model
    LOG << "Matching parameters from file with model:" << endl;
    vector<bool> usefile (ModelparamNames.size(),false);
    vector<int> oldloc (ModelparamNames.size(), 0);
    vector<bool> hasmatched (paramNames.size(), false); // to store if the file paramers have been matched
    for (unsigned p=0; p<ModelparamNames.size(); p++) {
      usefile[p]=false;
        for (unsigned q=0; q<paramNames.size(); q++) {
	  if (ModelparamNames[p] == paramNames[q]) {
	    usefile[p]=true;
	    oldloc[p] = q;
	    hasmatched[q]=true;
	    LOG << ModelparamNames[p] << ": Matched with file" << endl;
	  }
	}
	if (!usefile[p]) {
	  LOG << ModelparamNames[p] << ": Not matched, set from model default" << endl;
	}
    }

    //Make a note of any parameters in the file that were not matched
    for (unsigned int q=0; q<paramNames.size(); q++) {
      if (!hasmatched[q]) {
	LOG_ERR(paramNames[q] + ": Not matched!");
      }
    }

    // for (int a=0; a<usefile.size(); a++) {
//       cout << usefile[a] << "  " << oldloc[a] << endl;
//     }
 
    ColumnVector filemeans;
    ColumnVector modelmeans;
    ColumnVector newmeans;
    modelmeans = tempposterior.means;
    newmeans = modelmeans;

    SymmetricMatrix filecov;
    SymmetricMatrix modelcov;
    SymmetricMatrix newcov;
    modelcov = tempposterior.GetCovariance();
    newcov = modelcov;

    MVNDist fwddist;
    MVNDist noisedist;
    MVNDist newfwd(nmodparams);

    int nfwdparams = paramNames.size(); //number of fwd params in the MVN file
    int nnoiseparams = MVNfile[1]->means.Nrows() - nfwdparams; //number of noise parameters in the MVN file
    int nvox = MVNfile.size();

    for (int v=0; v<nvox; v++) {
      fwddist = MVNfile[v]->GetSubmatrix(1,nfwdparams);
      noisedist = MVNfile[v]->GetSubmatrix(nfwdparams+1,nfwdparams+nnoiseparams);

      for (unsigned int p=0; p<ModelparamNames.size(); p++) {
	// deal with the means
	if (usefile[p]) {
	  newmeans(p+1) = fwddist.means(oldloc[p]+1);
	}
      }
      newfwd.means=newmeans;

      //deal with the covariances
      filecov =fwddist.GetCovariance();
      for (unsigned int p=0; p<ModelparamNames.size(); p++) {
	for (int q=0; q<=p; q++) {
	  if(usefile[p]) {
	    if (usefile[q]) {
	      newcov(p+1,q+1) = filecov(oldloc[p]+1,oldloc[q]+1);
	    }
	  }
	}
      }
      newfwd.SetCovariance(newcov);

      //MVNDist distout(newfwd,noisedist);
      //continueFromDists.push_back(&distout);
      continueFromDists.push_back(new MVNDist(newfwd,noisedist) );
    }

  }

  

#endif //__FABBER_LIBRARYONLY
}

InferenceTechnique::~InferenceTechnique() 
{ 
  delete model;
  delete noise;
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

#include "inference_vb.h"
#include "inference_spatialvb.h"
#include "inference_nlls.h"

// Whenever you add a new class to inference.h, update this too.
InferenceTechnique* InferenceTechnique::NewFromName(const string& method)
{
  Tracer_Plus tr("PickInferenceTechnique");

  if (method == "vb")
    {
      return new VariationalBayesInferenceTechnique;
    }
  else if (method == "spatialvb")
    {
      return new SpatialVariationalBayes;
    }
  else if (method == "nlls")
    {
      return new NLLSInferenceTechnique;
    }
  else if (method == "")
    {
      throw Invalid_option("Must include the --method=vb or --method=spatialvb option");
    }
  else 
    {
      throw Invalid_option("Unrecognized --method: " + method);
    }
}


#ifdef __FABBER_MOTION

MCobj::MCobj(const DataSet& allData, int dof) {
  Tracer_Plus tr("MCobj::MCobj");

  //initialise
  mask = allData.GetMask();
  userdof=dof;
  num_iter=10;
  // the following sets up an initial zero deformation field
  Matrix datamat = allData.GetVoxelData();
  wholeimage.setmatrix(datamat,mask);
  modelpred=wholeimage;
  modelpred=0.0f;
  if (userdof > 12) {
    defx=modelpred;
    defx.setROIlimits(0,2);
    defx.activateROI();
    defx=defx.ROI();
    defy=defx;
    defz=defx;
    // Unnecessary initialisations?!?
    tmpx=defx;
    tmpy=defx;
    tmpz=defx;
  } else {
    mcf.setparams("verbose",false);
  }
  affmat=IdentityMatrix(4);
  finalimage=modelpred;
}


void MCobj::run_mc(const Matrix& modelpred_mat, Matrix& finalimage_mat) {
  Tracer_Plus tr("MCobj::run_mc");

  modelpred.setmatrix(modelpred_mat,mask);
  if (userdof>12) {
    UpdateDeformation(wholeimage,modelpred,num_iter,defx,defy,defz,finalimage,tmpx,tmpy,tmpz);
    defx=tmpx;
    defy=tmpy;
    defz=tmpz;
  } else {
    // mcf.register_volumes(4D reference,4D input image,refweight,inweight,4D output image);
    affmat = mcf.register_volumes(modelpred,wholeimage,mask,mask,finalimage);
    // apply transforms to wholeimage to get finalimage (the above is a dummy)
    for (int n=0; n<=wholeimage.maxt(); n++) {
      affine_transform(wholeimage[n],finalimage[n],affmat.SubMatrix(4*n+1,4*n+4,1,4));
    }
  }
  finalimage_mat = finalimage.matrix(mask);
}

#endif //__FABBER_MOTION
