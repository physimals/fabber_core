/*  fwdmodel_asl_2cpt.cc - Implements 2-compartment resting state ASL model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2009 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_2cpt.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
#include "miscmaths/miscprob.h"
using namespace NEWIMAGE;
#include "easylog.h"

string TwoCptFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_asl_grase.cc,v 1.14 2008/07/28 15:01:43 chappell Exp $";
}

void TwoCptFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("TwoCptFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // Tissue bolus perfusion
     prior.means(tiss_index()) = 0;
     precisions(tiss_index(),tiss_index()) = 1e-12;

     
     if (!singleti) {
       // Tissue bolus transit delay
       prior.means(tiss_index()+1) = 0.7;
       precisions(tiss_index()+1,tiss_index()+1) = 10;
     }
    
    
    // Tissue bolus length
     if (infertau) {
       prior.means(tau_index()) = seqtau;
       precisions(tau_index(),tau_index()) = 10;
     }

     if (infertaub)
       {
	 prior.means(taub_index()) = seqtau;
	 precisions(taub_index(),taub_index()) = 10;
       }

    // Arterial Perfusion & bolus delay

    if (inferart)
      {
	int aidx = art_index();
	prior.means(aidx) = 0;
	prior.means(aidx+1) = 0.5;
	precisions(aidx+1,aidx+1) = 10;
	precisions(aidx,aidx) = 1e-12;
      }
 
    // T1 & T1b
    if (infert1) {
      int tidx = t1_index();
      prior.means(tidx) = t1;  
      prior.means(tidx+1) = t1b; 
      precisions(tidx,tidx) = 100;
      precisions(tidx+1,tidx+1) = 100;
    }

    if (inferPS) {
      prior.means(PS_index()) = 30;
      precisions(PS_index(),PS_index()) = 0.01;
    }

    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital posterior
    // Tissue perfusion
    posterior.means(tiss_index()) = 10;
    precisions(tiss_index(),tiss_index()) = 1;
    // Arterial perfusion
    if (inferart)
      {
	posterior.means(art_index()) = 10;
	precisions(art_index(),art_index()) = 1;
      }
    posterior.SetPrecisions(precisions);
    
}    
    
    

void TwoCptFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("TwoCptFwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
   for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
      }
  
     // sensible limits on transit times
   if (!singleti) {
  if (params(tiss_index()+1)>timax-0.2) { paramcpy(tiss_index()+1) = timax-0.2; }
   }
  if (inferart) {
    if (params(art_index()+1)>timax-0.2) { paramcpy(art_index()+1) = timax-0.2; }
  }


  // parameters that are inferred - extract and give sensible names
  float ftiss;
  float delttiss;
  float tauset; //the value of tau set by the sequence (may be effectively infinite)
  float taubset;
  float fblood;
  float deltblood;
  float T_1;
  float T_1b;
  //  float inveffslope;
  //float trailingperiod;

  ftiss=paramcpy(tiss_index());
  if (!singleti) {
  delttiss=paramcpy(tiss_index()+1);
  }
  else {
    //only inferring on tissue perfusion, assume fixed value for tissue arrival time
    delttiss = 0;
  }

  if (infertau) { 
    tauset=paramcpy(tau_index()); 
  }
  else { tauset = seqtau;
  }
  
  if (infertaub) {
    taubset = paramcpy(taub_index());
  }
  else
    { taubset = tauset; }

  if (inferart) {
    fblood=paramcpy(art_index());
    deltblood=paramcpy(art_index()+1);
  }
  else {
    fblood = 0;
    deltblood = 0;
  }

  if (infert1) {
    T_1 = paramcpy(t1_index());
    T_1b = paramcpy(t1_index()+1);

    //T1 cannot be zero!
    if (T_1<1e-12) T_1=0.01;
    if (T_1b<1e-12) T_1b=0.01;
  }
  else {
    T_1 = t1;
    T_1b = t1b;
  }

  /*float bollen; //this is the bolus length as determined by the slope in inversion efficiency (onyl used with --inferinveff)
  if (inferinveff) { inveffslope = paramcpy(inveff_index());  bollen = delttiss+1/inveffslope; }
  else { inveffslope = 0; bollen = 1000; }  

  if (infertrailing) { trailingperiod = paramcpy(trailing_index()); }
  else { trailingperiod = 1e-6; }
  */


    float lambda = 0.9;

    float f_calib;
    // if we are using calibrated data then we can use ftiss to calculate T_1app
    if (calib) f_calib = ftiss;
    else       f_calib = 0.01; //otherwise assume sensible value (units of s^-1)

    float T_1app = 1/( 1/T_1 + f_calib/lambda );
    float R = 1/T_1app - 1/T_1b;

    float tau; //bolus length as seen by kintic curve
    float taub; //bolus length of blood as seen in signal
    

    float F=0;
    float kctissue;
    float kcblood;


    // loop over tis
    float ti;
    result.ReSize(tis.Nrows()*repeats);

    for(int it=1; it<=tis.Nrows(); it++)
      {
	ti = tis(it);
	F = 2*ftiss * exp(-ti/T_1app);

	/* According to EAGLE GRASE sequence bolus length is current TI - 0.1s (assuming infite length 'true' bolus)
	   However, also allow here bolus length to be finite as recorded in tauset
	   NB tauset is the 'true' bolus length, tau is what the tissue actually sees as a result of the sequence
	   25-3-2009 now deal with this scenario via pretisat parameter
	*/

	/*	if (grase)
	  {
	  //GRASE -  deal with bolus length (see above) */

	// Deal with saturation of the bolus before the TI - defined by pretisat
	if(tauset < ti - pretisat)
	  { tau = tauset; }
	else
	  { tau = ti -  pretisat; }
	
	if(taubset < ti -  pretisat)
	  {taub = taubset; }
	else
	  {taub = ti -  pretisat; }

	
	// tissue contribution
	// Using 2cpt model of Parkes and Tofts
	//Taking simplified solution - neglecting backflow
	
	if (!slow) { //fast model
	    // various variables
	    float PS=1;
	    float cbv=0.01; //NB not strictly cbv since we also have conc water in bloog to account for.
	    float kct_b; float kct_e;
	    float f=ftiss;

	    if(ti < delttiss)
		  { kctissue = 0;}
	    else {
	      float t;
	      if(ti >= delttiss && ti <= (delttiss + tau)) t = ti;
	      else t = delttiss+tau;

	      kct_b = 2.0*T_1b*cbv*f*(-exp(-(t*cbv-delttiss*T_1b*PS-delttiss*cbv-delttiss*T_1b*f_calib+t*T_1b*PS+t*cbv+t*T_1b*f_calib)/T_1b/cbv)+exp(-t/T_1b))/(T_1b*PS+cbv+T_1b*f_calib);

	      float MapleGenVar1 = -2.0*T_1b*f*T_1;      
	      float MapleGenVar4 = PS;      
	      float MapleGenVar6 = -PS*T_1b*T_1*exp(-(t*T_1-t*T_1b)/T_1b/T_1)+PS*T_1b*T_1*exp((-t*T_1+delttiss*T_1b)/T_1b/T_1)+T_1b*T_1*exp((-t*T_1+delttiss*T_1b)/T_1b/T_1)*f_calib-T_1b*T_1*exp(-(t*T_1-t*T_1b)/T_1b/T_1)*f_calib-T_1b*cbv*exp(-(t*cbv*T_1+t*PS*T_1b*T_1+t*cbv*T_1+t*f_calib*T_1b*T_1-t*cbv*T_1b-delttiss*T_1*T_1b*PS-delttiss*T_1*cbv-delttiss*T_1*T_1b*f_calib)/cbv/T_1b/T_1);      
	      float MapleGenVar5 = MapleGenVar6+T_1b*exp(-(t*T_1-t*T_1b)/T_1b/T_1)*cbv+T_1b*exp(-(t*cbv*T_1-delttiss*T_1*T_1b*PS-delttiss*T_1*cbv-delttiss*T_1*T_1b*f_calib+delttiss*PS*T_1b*T_1+delttiss*T_1*cbv+delttiss*f_calib*T_1b*T_1-delttiss*cbv*T_1b)/T_1b/T_1/cbv)*cbv-T_1b*exp((-t*T_1+delttiss*T_1b)/T_1b/T_1)*cbv-exp(-(t*T_1-t*T_1b)/T_1b/T_1)*T_1*cbv+exp((-t*T_1+delttiss*T_1b)/T_1b/T_1)*cbv*T_1;      
	      float MapleGenVar3 = MapleGenVar4*MapleGenVar5;      
	      MapleGenVar4 = exp(-1/T_1*t)/(T_1b*T_1b*PS*PS*T_1+2.0*T_1b*PS*T_1*cbv+2.0*T_1b*T_1b*PS*f_calib*T_1-T_1b*T_1b*PS*cbv+T_1*cbv*cbv+2.0*cbv*f_calib*T_1b*T_1-T_1b*cbv*cbv+T_1b*T_1b*f_calib*f_calib*T_1-T_1b*T_1b*f_calib*cbv);      
	      float MapleGenVar2 = MapleGenVar3*MapleGenVar4;
	      kct_e = MapleGenVar1*MapleGenVar2;

	      if(ti >= delttiss && ti <= (delttiss + tau)) kctissue = kct_e+kct_b;
	      else {
		float kct2_b = kct_b*exp((T_1b*PS+cbv+T_1b*f_calib)*(delttiss+tau-ti)/T_1b/cbv);
		float X = kct_b; float Y=kct_e;
		float kct2_e = (-PS*T_1b*T_1*X*exp((-ti*PS*T_1b*T_1-ti*cbv*T_1-ti*f_calib*T_1b*T_1+ti*cbv*T_1b+delttiss*PS*T_1b*T_1+T_1*T_1b*PS*tau+delttiss*T_1*cbv+T_1*cbv*tau+delttiss*f_calib*T_1b*T_1+T_1*T_1b*f_calib*tau)/cbv/T_1b/T_1)+exp(1/T_1*(delttiss+tau))*PS*T_1b*T_1*X+exp(1/T_1*(delttiss+tau))*Y*PS*T_1b*T_1+exp(1/T_1*(delttiss+tau))*Y*T_1*cbv+exp(1/T_1*(delttiss+tau))*Y*f_calib*T_1b*T_1-exp(1/T_1*(delttiss+tau))*Y*T_1b*cbv)*exp(-1/T_1*ti)/(PS*T_1b*T_1+T_1*cbv+f_calib*T_1b*T_1-T_1b*cbv);

		kctissue = kct2_e + kct2_b;
	      }
	    }
	}
	else  { // slow model
    
	  float PScbv;
	  // Deal with PS - NB define PScbv = PS/cbv
	  if (inferPS) {
	    PScbv=paramcpy(PS_index());
	  }
	  else PScbv = 30;
    
	  float kct_b; float kct_e;
	    
	  if(ti < delttiss)
	    { kctissue = 0;}
	  else {
	    float t;
	    if(ti >= delttiss && ti <= (delttiss + tau)) t = ti;
	    else t = delttiss+tau;

	    kct_b = 2.0*ftiss*(exp(-t/T_1b)-exp((-t-t*PScbv*T_1b+delttiss*PScbv*T_1b)/T_1b))/PScbv;
	    kct_e = -2.0*T_1*T_1b*ftiss*(T_1b*T_1*exp(delttiss*(T_1b-T_1)/T_1/T_1b)*PScbv-PScbv*T_1b*exp(t*(T_1b-T_1)/T_1/T_1b)*T_1-T_1b*exp((-t*T_1-t*PScbv*T_1b*T_1+t*T_1b+delttiss*PScbv*T_1b*T_1)/T_1b/T_1)+T_1b*exp(t*(T_1b-T_1)/T_1/T_1b)-exp(t*(T_1b-T_1)/T_1/T_1b)*T_1+exp((-t*T_1-t*PScbv*T_1b*T_1+t*T_1b+delttiss*PScbv*T_1b*T_1)/T_1b/T_1)*T_1)*exp(-1/T_1*t)/(2.0*T_1b*T_1+PScbv*T_1b*T_1b*T_1-T_1b*T_1b-T_1*T_1-PScbv*T_1b*T_1*T_1);
	
	    if(ti >= delttiss && ti <= (delttiss + tau)) kctissue = kct_e+kct_b;
	    else {
	      float X = kct_b; float Y=kct_e;
	      float kct2_b = X*exp((1.0+PScbv*T_1b)*(-t+delttiss+tau)/T_1b);
	      float kct2_e = (-PScbv*T_1b*T_1*X*exp((-t*T_1-t*PScbv*T_1b*T_1+t*T_1b+delttiss*T_1+T_1*tau+delttiss*PScbv*T_1b*T_1+T_1*PScbv*T_1b*tau)/T_1b/T_1)+exp(1/T_1*(delttiss+tau))*PScbv*T_1b*T_1*X+exp(1/T_1*(delttiss+tau))*Y*T_1+exp(1/T_1*(delttiss+tau))*Y*PScbv*T_1b*T_1-exp(1/T_1*(delttiss+tau))*Y*T_1b)*exp(-1/T_1*t)/(T_1+PScbv*T_1b*T_1-T_1b);
	kctissue = kct2_e + kct2_b;
	      }
	    }
	}


	if (kctissue<0) { kctissue = 0; } //dont allow negative values (should be redundant)

	    /* arterial contribution */
	    if(ti < deltblood)
	      { 
		//kcblood = 0;
		/* use a artifical lead in period for arterial bolus to improve model fitting */
		kcblood = fblood * exp(-deltblood/T_1b) * (0.98 * exp( (ti-deltblood)/0.05 ) + 0.02 * ti/deltblood );
	      }
	    else if(ti >= deltblood && ti <= (deltblood + taub))
	      { 
		kcblood = fblood * exp(-ti/T_1b); 
		if (kcblood<0) { kcblood = 0; } //dont allow negative values
	      }
	    else /*(ti > deltblood + tau) */
	      {
		kcblood = 0; //end of bolus
		/* artifical lead out period for taub model fitting */
		kcblood = fblood * exp(-(deltblood+taub)/T_1b)  * (0.98 * exp( -(ti - deltblood - taub)/0.05) + 0.02 * (1-(ti - deltblood - taub)/5));
		if (kcblood<0) kcblood=0; //negative values are possible with the lead out period equation
									   
	      }

	    if (isnan(kctissue)) { kctissue=0; LOG << "Warning NaN in tissue curve at TI:" << ti << " with f:" << ftiss << " delt:" << delttiss << " tau:" << tau << " T1:" << T_1 << " T1b:" << T_1b << endl; }
	    //}

	/* output */
	// loop over the repeats
	for (int rpt=1; rpt<=repeats; rpt++)
	  {
	    result( (it-1)*repeats+rpt ) = kctissue + kcblood;
	  }

 
      }
    //cout << result.t();
    

  return;
}

TwoCptFwdModel::TwoCptFwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
      repeats = convertTo<int>(args.Read("repeats")); // number of repeats in data
      t1 = convertTo<double>(args.ReadWithDefault("t1","1.3"));
      t1b = convertTo<double>(args.ReadWithDefault("t1b","1.5"));

      pretisat = convertTo<double>(args.ReadWithDefault("pretisat","0")); // deal with saturation of the bolus a fixed time pre TI measurement
      grase = args.ReadBool("grase"); // DEPRECEATED data has come from the GRASE-ASL sequence - therefore apply pretisat of 0.1s
      if (grase) pretisat=0.1;

      calib = args.ReadBool("calib");
      slow = !(args.ReadBool("fast")); //option to use 'fast solution'
      inferPS = true;

      infertau = args.ReadBool("infertau"); // infer on bolus length?
      infert1 = args.ReadBool("infert1"); //infer on T1 values?
      inferart = args.ReadBool("inferart"); //infer on arterial compartment?
      //inferinveff = args.ReadBool("inferinveff"); //infer on a linear decrease in inversion efficiency?
      //infertrailing = args.ReadBool("infertrailing"); //infers a trailing edge bolus slope using new model
      seqtau = convertTo<double>(args.ReadWithDefault("tau","1000")); //bolus length as set by sequence (default of 1000 is effectively infinite
      bool ardoff = false;
      ardoff = args.ReadBool("ardoff");
      bool tauboff = false;
      tauboff = args.ReadBool("tauboff"); //forces the inference of arterial bolus off

      // combination options
      infertaub = false;
      if (inferart && infertau && !tauboff) infertaub = true;

      // deal with ARD selection
      doard=false;
      if (inferart==true && ardoff==false) { doard=true; }

      
      /* if (infertrailing) {
	if (!infertau) {
	  // do not permit trailing edge inference without inferring on bolus length
	  throw Invalid_option("--infertrailing has been set without setting --infertau");
	}
	else if (inferinveff)
	  //do not permit trailing edge inference and inversion efficiency inference (they are mututally exclusive)
	  throw Invalid_option("--infertrailing and --inferinveff may not both be set");
	  }*/

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
      
      singleti = false; //normally we do multi TI ASL
      if (tis.Nrows()==1) {
	//only one TI therefore only infer on CBF and ignore other inference options
	LOG << "--Single inversion time mode--" << endl;
	LOG << "Only a sinlge inversion time has been supplied," << endl;
	LOG << "Therefore only tissue perfusion will be inferred." << endl;
	LOG << "-----" << endl;
	singleti = true;
	// force other inference options to be false
	infertau = false; infert1 = false; inferart = false; //inferinveff = false;
      }
	
      // add information about the parameters to the log
      LOG << "Inference using Buxton Kinetic Curve model" << endl;
      if (pretisat>0) LOG << "Saturation of" << pretisat << "s before TI has been specified" << endl;
      if (grase) LOG << "Using pre TI saturation of 0.1 for GRASE-ASL sequence" << endl;
      if (calib) LOG << "Input data is in physioligcal units, using estimated CBF in T_1app calculation" << endl;
      LOG << "    Data parameters: #repeats = " << repeats << ", t1 = " << t1 << ", t1b = " << t1b;
      LOG << ", bolus length (tau) = " << seqtau << endl ;
      if (infertau) {
	LOG << "Infering on bolus length " << endl; }
      if (inferart) {
	LOG << "Infering on artertial compartment " << endl; }
      if (doard) {
	LOG << "ARD has been set on arterial compartment " << endl; }
      if (infert1) {
	LOG << "Infering on T1 values " << endl; }
      /*if (inferinveff) {
	LOG << "Infering on Inversion Efficency slope " << endl; }
      if (infertrailing) {
      LOG << "Infering bolus trailing edge period" << endl; }*/
      LOG << "TIs: ";
      for (int i=1; i <= tis.Nrows(); i++)
	LOG << tis(i) << " ";
      LOG << endl;
	  
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    
    
 
}

void TwoCptFwdModel::ModelUsage()
{ 
  cout << "\nUsage info for --model=grase:\n"
       << "Required parameters:\n"
       << "--repeats=<no. repeats in data>\n"
       << "--ti1=<first_inversion_time_in_seconds>\n"
       << "--ti2=<second_inversion_time>, etc...\n"
       << "Optional arguments:\n"
       << "--grase *DEPRECEATAED* (data collected using GRASE-ASL: same as --pretissat=0.1)\n"
       << "--pretisat=<presat_time> (Define that blood is saturated a specific time before TI image acquired)\n"
       << "--calib (data has been provided in calibrated units)\n"
       << "--tau=<temporal_bolus_length> (default 10s if --infertau not set)\n"
       << "--t1=<T1_of_tissue> (default 1.3)\n"
       << "--t1b=<T1_of_blood> (default 1.5)\n"
       << "--infertau (to infer on bolus length)\n"
       << "--inferart (to infer on arterial compartment)\n"
       << "--infert1 (to infer on T1 values)\n"
    ;
}

void TwoCptFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void TwoCptFwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  names.push_back("ftiss");
  if (!singleti) 
    names.push_back("delttiss");
  if (infertau)
    {
    names.push_back("tautiss");
    }
  if (inferart) {
    names.push_back("fblood");
    names.push_back("deltblood");
  }
  if (infert1) {
    names.push_back("T_1");
    names.push_back("T_1b");
  }
  /* if (inferinveff) {
    names.push_back("Inveffslope");
  }
  if (infertrailing) {
    names.push_back("trailingperiod");
    }*/
  if (infertaub) {
    names.push_back("taublood");
  }
  if (inferPS) {
    names.push_back("PS_cbv");
  }
}

void TwoCptFwdModel::SetupARD( const MVNDist& theta, MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("TwoCptFwdModel::SetupARD");

 

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

void TwoCptFwdModel::UpdateARD(
				const MVNDist& theta,
				MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("TwoCptFwdModel::UpdateARD");
  
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
