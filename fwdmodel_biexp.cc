/*  fwdmodel_biexp.cc - Implements a simple bi-exponential model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_biexp.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

string BiExpFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_biexp.cc,v 1.13 2008/04/10 15:06:27 chappell Exp $";
}

void BiExpFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("BiExpFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors

    // 1st Exp decay
     prior.means(1) = 0;
     precisions(1,1) = 1e-12;
     prior.means(2) = -1; //this is exp
     precisions(2,2) = 1;
    
     // 2nd Exp decay
     prior.means(3) = 0;
     precisions(3,3) = 1e-12;
     prior.means(4) = 1; //this is exp
     precisions(4,4) = 1;

    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;

    
}    
    
    

void BiExpFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("BiExpFwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
   for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
      }

  // parameters that are inferred - extract and give sensible names
   float mag1;
   float tc1;
   float mag2;
   float tc2;

   mag1=params(1);
   tc1=exp(params(2));
   mag2=params(3);
   tc2=exp(params(4));

   float imax = 51;
   float t0 = 0;
   float tint = 0.1;

   result.ReSize(imax);
   // loop to create result
   for (int i=1; i<=imax; i++)
     {
       float ti = t0 + (i-1)*tint;
       result(i) = mag1*exp(-tc1*ti) + mag2*exp(-tc2*ti);
 
      }
    //cout << result.t();
    

  return;
}

BiExpFwdModel::BiExpFwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
      //repeats = convertTo<int>(args.Read("repeats")); // number of repeats in data
      

      // deal with ARD selection
      //doard=false;
      //if (inferart==true && ardoff==false) { doard=true; }

     
      // add information about the parameters to the log
      LOG << "Inference using Bi-exponetial model" << endl;
      
	  
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    
    
 
}

void BiExpFwdModel::ModelUsage()
{ 
  cout << "\nUsage info for --model=biexp:\n"
       << "Required parameters:\n"
       
       << "Optional arguments:\n"
       
    ;
}

void BiExpFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void BiExpFwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  names.push_back("mag1");
  names.push_back("tc1");
  names.push_back("mag2");
  names.push_back("tc2");
}

void BiExpFwdModel::SetupARD( const MVNDist& theta, MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("BiExpFwdModel::SetupARD");

  int ardindex = ard_index();


   if (doard)
    {
      SymmetricMatrix PriorPrec;
      PriorPrec = thetaPrior.GetPrecisions();
      
      PriorPrec(ardindex,ardindex) = 1e-12; //set prior to be initally non-informative
      
      thetaPrior.SetPrecisions(PriorPrec);

      thetaPrior.means(ardindex)=0;

      //set the Free energy contribution from ARD term
      SymmetricMatrix PostCov = theta.GetCovariance();
      double b = 2/(theta.means(ardindex)*theta.means(ardindex) + PostCov(ardindex,ardindex));
      Fard = -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
    }

  return;
}

void BiExpFwdModel::UpdateARD(
				const MVNDist& theta,
				MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("BiExpFwdModel::UpdateARD");
  
  int ardindex = ard_index();


  if (doard)
    {
      SymmetricMatrix PriorCov;
      SymmetricMatrix PostCov;
      PriorCov = thetaPrior.GetCovariance();
      PostCov = theta.GetCovariance();

      PriorCov(ardindex,ardindex) = theta.means(ardindex)*theta.means(ardindex) + PostCov(ardindex,ardindex);

      
      thetaPrior.SetCovariance(PriorCov);

      //Calculate the extra terms for the free energy
      double b = 2/(theta.means(ardindex)*theta.means(ardindex) + PostCov(ardindex,ardindex));
      Fard = -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
    }

  return;

  }
