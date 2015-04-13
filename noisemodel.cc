/*  noisemodel.cc - Class implementation for generic noise models

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "noisemodel.h"

// Handy mathematical function, used in some free energy calculations
#include <math.h>

// Calculate log-gamma from a Taylor expansion; good to one part in 2e-10.
double gammaln(double x)
{
  ColumnVector series(7);
  series << 2.5066282746310005 << 76.18009172947146
	 << -86.50532032941677 << 24.01409824083091
	 << -1.231739572450155 << 0.1208650973866179e-2
	 << -0.5395239384953e-5;

  double total = 1.000000000190015;
  for (int i = 2; i <= series.Nrows(); i++)
    total += series(i) / (x+i-1);

  return log( series(1) * total / x ) + (x + 0.5)*log(x + 5.5) - x - 5.5;
}

#include "noisemodel_ar.h"
#include "noisemodel_white.h"

NoiseModel* NoiseModel::NewFromName(const string& name, ArgsType& args)
{
    // Update this to add your own models to the code
    
  if (name == "ar")
    {
      string nPhis = args.ReadWithDefault("num-echoes","(default)");
      if (nPhis == "(default)")
	{
	  nPhis = "1";
	  Warning::IssueOnce("Defaulting to --num-echoes=1");
	}
      string ar1CrossTerms = args.ReadWithDefault("ar1-cross-terms","none");
      
      return new Ar1cNoiseModel(ar1CrossTerms, convertTo<int>(nPhis));
    }
  else if (name == "ar1" || name == "ar1c")
    {
      Warning::IssueOnce("--noise="+name+" is depreciated; use --noise=ar --num-echoes=2 for dual-echo ASL.");
      string ar1CrossTerms = args.ReadWithDefault("ar1-cross-terms","none");
      string nPhis = args.ReadWithDefault("num-echoes","2");
      
      return new Ar1cNoiseModel(ar1CrossTerms, convertTo<int>(nPhis));
    }
    else if (name == "white")
    {
      //      string pattern = args.ReadWithDefault("noise-pattern","1");
      //      return new WhiteNoiseModel(pattern);
      return new WhiteNoiseModel(args);
    }
    // Your models go here!
    else
    {
      throw Invalid_option("Unrecognized noise model '" + name + "'");
    }
}

// ARD stuff
double NoiseModel::SetupARD(vector<int> ardindices,
			  const MVNDist& theta,
			  MVNDist& thetaPrior) const {
  Tracer_Plus tr("Noisemodel::SetupARD");
  double Fard=0;

  if (~ardindices.empty()) {
    SymmetricMatrix PriorPrec;
    PriorPrec = thetaPrior.GetPrecisions();
    SymmetricMatrix PostCov = theta.GetCovariance();

    for (int i=0; i< ardindices.size(); i++)
      {
	PriorPrec(ardindices[i],ardindices[i]) = 1e-12; //set prior to be initally non-informative
	thetaPrior.means(ardindices[i]) = 0;
	
	//set the Free energy contribution from ARD term
	double b = 2/(theta.means(ardindices[i])*theta.means(ardindices[i]) + PostCov(ardindices[i],ardindices[i]));
	Fard += -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
      }

	thetaPrior.SetPrecisions(PriorPrec);
  }

  return Fard;
}

double NoiseModel::UpdateARD(vector<int> ardindices,
			  const MVNDist& theta,
			  MVNDist& thetaPrior) const {
  Tracer_Plus tr("Noisemodel::UpdateARD");
  double Fard=0;

  if (~ardindices.empty()) {
    SymmetricMatrix PriorCov;
      SymmetricMatrix PostCov;
      PriorCov = thetaPrior.GetCovariance();
      PostCov = theta.GetCovariance();

    for (int i=0; i< ardindices.size(); i++)
      {
	PriorCov(ardindices[i],ardindices[i]) = theta.means(ardindices[i])*theta.means(ardindices[i]) + PostCov(ardindices[i],ardindices[i]);
	
	//set the Free energy contribution from ARD term
	double b = 2/(theta.means(ardindices[i])*theta.means(ardindices[i]) + PostCov(ardindices[i],ardindices[i]));
	Fard += -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
      }

	thetaPrior.SetCovariance(PriorCov);
  }

  return Fard;
}
