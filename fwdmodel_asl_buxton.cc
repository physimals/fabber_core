/*  fwdmodel_asl_buxton.cc - Implements the Buxton kinetic curve model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_buxton.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

string BuxtonFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_asl_buxton.cc,v 1.8 2010/02/11 16:48:42 chappell Exp $";
}

void BuxtonFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("BuxtonFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // Tissue bolus perfusion
     prior.means(tiss_index()) = 0;
     precisions(tiss_index(),tiss_index()) = 1e-12;

    // Tissue bolus transit delay
     prior.means(tiss_index()+1) = 0.7;
     precisions(tiss_index()+1,tiss_index()+1) = 10;
    
    
    // Tissue bolus length
     if (infertau) {
       prior.means(tau_index()) = seqtau;
       precisions(tau_index(),tau_index()) = 100;
     }

     // T1 & T1b
    if (infert1) {
      int tidx = t1_index();
      prior.means(tidx) = t1;  
      prior.means(tidx+1) = t1b; 
      precisions(tidx,tidx) = 100;
      precisions(tidx+1,tidx+1) = 100;
    }

    // second bolus
    if (twobol) {
      int tiss2idx = tiss2_index();
      // perfusion
      prior.means(tiss2idx) = 0;
      precisions(tiss2idx,tiss2idx) = 1e-12; //should be largely irrelevant as this will be subject to ARD

      // transit delay
      prior.means(tiss2idx+1) = 1;
      precisions(tiss2idx+1,tiss2idx+1) = 10;
    }

    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;

    // Modify posterior for perfusion to a more realisitic starting point
    posterior.means(tiss_index()) = 10;
    precisions(tiss_index(),tiss_index()) = 0.1;
    posterior.SetPrecisions(precisions);
    
}    
    
    

void BuxtonFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("BuxtonFwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
  for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
    }

     // sensible limits on transit times
     if (params(tiss_index()+1)>timax-0.2) { paramcpy(tiss_index()+1) = timax-0.2; }
 
  // parameters that are inferred - extract and give sensible names
  float ftiss;
  float delttiss;
  float tautiss;
  float T_1;
  float T_1b;

  float ftiss2 = 0;
  float delttiss2 = 0;

  ftiss=paramcpy(tiss_index());
  delttiss=paramcpy(tiss_index()+1);

  if (infertau) { tautiss=paramcpy(tau_index()); }
  else { tautiss = seqtau; }
  
  if (infert1) {
    T_1 = paramcpy(t1_index());
    T_1b = paramcpy(t1_index()+1);
  }
  else {
    T_1 = t1;
    T_1b = t1b;
  }

  if (twobol) {
    ftiss2 = paramcpy(tiss2_index());
    delttiss2 = paramcpy(tiss2_index()+1);
  }
    

  //float lambda = 0.9;

    float T_1app = 1/( 1/T_1 + 0.01/lambda );
    float R = 1/T_1app - 1/T_1b;

    float tau1 = delttiss;
    float tau2 = delttiss + tautiss;

    // for second bolus
    float tau3 = delttiss2;
    float tau4 = delttiss2 + tautiss;

    float F=0;float F2=0;
    float kctissue; float kctissue2;
 

    // loop over tis
    float ti;
    result.ReSize(tis.Nrows()*repeats);

    for(int it=1; it<=tis.Nrows(); it++)
      {
	ti = tis(it);

	// 1st bolus
	F = 2*ftiss * exp(-ti/T_1app);

	/* tissue contribution */
	if(ti < tau1)
	  { kctissue = 0;}
	else if(ti >= tau1 && ti <= tau2)
	  {
	    kctissue = F/R * (exp(R*ti) - exp(R*tau1));
	  }
	else /*(ti > tau2)*/
	  {kctissue = F/R * (exp(R*tau2) - exp(R*tau1)); }

	// 2nd bolus
	F2 = 2*ftiss2/lambda * exp(-ti/T_1app);

	/* tissue contribution */
	if(ti < tau3)
	  { kctissue2 = 0;}
	else if(ti >= tau3 && ti <= tau4)
	  {
	    kctissue2 = F2/R * (exp(R*ti) - exp(R*tau3));
	  }
	else /*(ti > tau2)*/
	  {kctissue2 = F2/R * (exp(R*tau4) - exp(R*tau3)); }
	
	/* output */
	// loop over the repeats
	for (int rpt=1; rpt<=repeats; rpt++)
	  {
	    result( (it-1)*repeats+rpt ) = kctissue + kctissue2;
	  }


      }


  return;
}

BuxtonFwdModel::BuxtonFwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
      repeats = convertTo<int>(args.Read("repeats")); // number of repeats in data
      t1 = convertTo<double>(args.ReadWithDefault("t1","1.3"));
      t1b = convertTo<double>(args.ReadWithDefault("t1b","1.5"));
      lambda = convertTo<double>(args.ReadWithDefault("lambda","0.9"));
      infertau = args.ReadBool("infertau"); // infer on bolus length?
      infert1 = args.ReadBool("infert1"); //infer on T1 values?
      twobol = args.ReadBool("twobol"); //infer a second bolus?
      doard=false;
      if (twobol) { doard = true; } //ARD is always used on perfusion of second bolus
      seqtau = convertTo<double>(args.Read("tau")); 
 
      // Deal with tis
      tis.ReSize(1); //will add extra values onto end as needed
      tis(1) = atof(args.Read("ti1").c_str());
      
      while (true) //get the rest of the tis
	{
	  int N = tis.Nrows()+1;
	  string tiString = args.ReadWithDefault("ti"+stringify(N), "stop!");
	  if (tiString == "stop!") break; //we have run out of tis
	 
	  // append the new ti onto the end of the list
	  ColumnVector tmp(1);
	  tmp = convertTo<double>(tiString);
	  tis &= tmp; //vertical concatenation

	}
      timax = tis.Maximum(); //dtermine the final TI

      // add information about the parameters to the log
      LOG << "    Data parameters: #repeats = " << repeats << ", t1 = " << t1 << ", t1b = " << t1b;
      LOG << ", bolus length (tau) = " << seqtau << endl ;
      if (infertau) {
	LOG << "Infering on bolus length " << endl; }
      if (infert1) {
	LOG << "Infering on T1 values " << endl; }
      if (twobol) {
	LOG << "Inferring upon a second bolus " << endl; }
      LOG << "TIs: ";
      for (int i=1; i <= tis.Nrows(); i++)
	LOG << tis(i) << " ";
      LOG << endl;
	  
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    
    
 
}

void BuxtonFwdModel::ModelUsage()
{ 
  cout << "\nUsage info for --model=buxton:\n"
       << "Required parameters:\n"
       << "--repeats=<no. repeats in data>\n"
       << "--ti1=<first_inversion_time_in_seconds>\n"
       << "--ti2=<second_inversion_time>, etc...\n"
       << "--tau=<temporal_bolus_length_in_seconds> \n"
       << "Optional arguments:\n"
       << "--t1=<T1_of_tissue> (default 1.3)\n"
       << "--t1b=<T1_of_blood> (default 1.5)\n"
       << "--infertau (to infer on bolus length)\n"
       << "--infert1 (to infer on T1 values)\n"
    ;
}

void BuxtonFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void BuxtonFwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  names.push_back("ftiss");
  names.push_back("delttiss");
  if (infertau)
    names.push_back("tautiss");
   if (infert1) {
    names.push_back("T_1");
    names.push_back("T_1b");
  }
   if (twobol) {
     names.push_back("ftiss2");
     names.push_back("delttiss2");
   }
}

/* ARD for second bolus perfusion */
/* taken from fwdmodel_asl_grase.cc (29-11-2007) */
void BuxtonFwdModel::SetupARD( const MVNDist& theta, MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("BuxtonFwdModel::SetupARD");

  int ardindex = ard_index();

   if (doard)
    {
      SymmetricMatrix PriorPrec;
      PriorPrec = thetaPrior.GetPrecisions();
      
      PriorPrec(ardindex,ardindex) = 1e-12;
      
      thetaPrior.SetPrecisions(PriorPrec);

      thetaPrior.means(ardindex)=0;

      //set the Free energy contribution from ARD term
      SymmetricMatrix PostCov = theta.GetCovariance();
      double b = 2/(theta.means(ardindex)*theta.means(ardindex) + PostCov(ardindex,ardindex));
      Fard = -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
    }

  return;
}

void BuxtonFwdModel::UpdateARD(
				const MVNDist& theta,
				MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("BuxtonFwdModel::UpdateARD");
  
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
