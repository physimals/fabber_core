/*  fwdmodel_biexp.cc - Implements a model for correcting off resonance effect for multiphase pcASL

    Michael Chappell, QuBIc (IBME) & FMRIB Image Analysis Group

    Copyright (C) 2013 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_multiphase.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

string MultiPhaseASLFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_asl_multiphase.cc,v 1.1 2013/08/08 12:30:47 chappell Exp $";
}

void MultiPhaseASLFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("MultiPhaseASLFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors

    // magnitude
     prior.means(1) = 0;
     precisions(1,1) = 1e-12;

     // phase (radians)
     prior.means(2) = 0;
     precisions(2,2) = M_PI/10;
    
     // offset
     prior.means(3) = 0;
     precisions(3,3) = 1e-12;

    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;

    
}    

void MultiPhaseASLFwdModel::Initialise(MVNDist& posterior) const
{
  Tracer_Plus tr("MultiPhaseASLFwdModel::Initialise");
  // init the magntidue and offset parameters

  // mean over the repeats
  ColumnVector dmean(8);
  dmean=0.0;
  for (int i=1; i<=8; i++) {
    for (int j=1; j<=repeats; j++) {
      dmean(i) = dmean(i) + data((j-1)*8+i);
    }
  }
  dmean = dmean/repeats;

  double dmax = dmean.Maximum();
  double dmin = dmean.Minimum();

  posterior.means(1) = (dmax-dmin)/2;
  posterior.means(3) = (dmax+dmin)/2;

//init the mid phase value value - by finding the point where the max intensity is
  int ind;
  float val;
  val = dmean.Maximum1(ind); // find the max
  val = (ind-1)*180/M_PI; //frequency of the minimum in ppm
  if (val>179) val -= 360;
  val *= M_PI/180;
  posterior.means(2) = val;
}
    
    

void MultiPhaseASLFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("MultiPhaseASLFwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
   for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
      }

  // parameters that are inferred - extract and give sensible names
   float mag;
   float phase;
   float offset;

   mag=params(1);
   phase=params(2)*180/M_PI;
   offset=params(3);

   int nn=8*repeats;
   result.ReSize(nn);
   // loop to create result
   for (int i=1; i<=8; i++)
     {
       double ph = (360/8*(i-1) );
       if (ph>179) ph -= 360;

       double evalfunc = mag*( -2/(1+exp( (abs(ph - phase) -55)/12)) ) + offset; //note using the given values requires phases here to be in degrees

       for (int j=1; j<=repeats; j++) {
	 result((j-1)*8+i) = evalfunc;
       }
 
      }
    //cout << result.t();

  return;
}

MultiPhaseASLFwdModel::MultiPhaseASLFwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
      repeats = convertTo<int>(args.Read("repeats")); // number of repeats in data

      // deal with ARD selection
      //doard=false;
      //if (inferart==true && ardoff==false) { doard=true; }

     
      // add information about the parameters to the log
      LOG << "Inference using Fermi model" << endl;
      
	  
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    
    
 
}

void MultiPhaseASLFwdModel::ModelUsage()
{ 
  cout << "\nUsage info for --model=biexp:\n"
       << "Required parameters:\n"
       
       << "Optional arguments:\n"
       
    ;
}

void MultiPhaseASLFwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  names.push_back("mag");
  names.push_back("phase");
  names.push_back("offset");
}
