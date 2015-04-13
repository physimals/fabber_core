/*  fwdmodel_asl_devel.cc - Developement resting state ASL model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2009 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_devel.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
#include "miscmaths/miscprob.h"
using namespace NEWIMAGE;
#include "easylog.h"

string DevelFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_asl_devel.cc,v 1.2 2010/10/01 13:46:54 chappell Exp $";
}

void DevelFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("DevelFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // Tissue bolus perfusion
     if (infertiss) {
     prior.means(tiss_index()) = 0;
     precisions(tiss_index(),tiss_index()) = 1e-12;

     
     //if (!singleti) {
       // Tissue bolus transit delay
       prior.means(tiss_index()+1) = 0.7;
       precisions(tiss_index()+1,tiss_index()+1) = 10;
       // }
    
     }

    // Tissue bolus length
     if (infertau && infertiss) {
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

    /* if (inferart) {
      prior.means(R_index()) = log(10);
      precisions(R_index(),R_index()) = 1;
      }*/

    if (inferwm) {
      int wmi = wm_index();
      prior.means(wmi) = 0;
      prior.means(wmi+1) = 1.2;
      precisions(wmi,wmi) = 1e-12;
      precisions(wmi+1,wmi+1) = 10;

      if (infertau) {
	prior.means(wmi+2) = seqtau;
	precisions(wmi+2,wmi+2) = 10;
      }

      if (infert1) {
	prior.means(wmi+3) = t1wm;
	precisions(wmi+3,wmi+3) = 100;
      }

      if (usepve) {
      //PV entries, the means get overwritten elsewhere if the right sort of prior is specified
      // default is to allow both (NB artifically defies sum(pve)=1)
      int pvi= pv_index();
      prior.means(pvi) = 1; //GM is first
      prior.means(pvi+1)= 1; //WM is second

      // precisions are big as we treat PV parameters as correct
      // NB they are not accesible from the data anyway
      precisions(pvi,pvi) = 1e12;
      precisions(pvi+1,pvi+1) = 1e12;
      }

    }

    //dispersion parameters
    prior.means(disp_index()) = 0.05;
    prior.means(disp_index()+1) = 3;
    precisions(disp_index(),disp_index()) = 400;
    precisions(disp_index()+1,disp_index()+1) = 1;

    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital posterior
    // Tissue perfusion
    if (infertiss) {
    posterior.means(tiss_index()) = 10;
    precisions(tiss_index(),tiss_index()) = 1;
    }
    // Arterial perfusion
    if (inferart)
      {
	posterior.means(art_index()) = 10;
	precisions(art_index(),art_index()) = 1;
      }

    if (inferwm)
      {
	posterior.means(wm_index()) = 10;
	precisions(wm_index(),wm_index()) = 1;
      }

    posterior.means(disp_index()) = 0.05;
    
    posterior.SetPrecisions(precisions);
    
}    
    
    

void DevelFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("DevelFwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
   for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
      }
  
     // sensible limits on transit times
   if (infertiss) {
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

  float pv_gm;
  float pv_wm;

  float fwm;
  float deltwm;
  float tauwmset;
  float T_1wm;

  float p;
  float s;

  //  float RR;
  //  float inveffslope;
  //float trailingperiod;

  if (infertiss) {
  ftiss=paramcpy(tiss_index());
  //if (!singleti) {
  delttiss=paramcpy(tiss_index()+1);
  //}
  //else {
    //only inferring on tissue perfusion, assume fixed value for tissue arrival time
    //delttiss = 0;
    //}
  }
  else {
    ftiss=0;
    delttiss=0;
  }

  if (infertau && infertiss) { 
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
    if (T_1<0.01) T_1=0.01;
    if (T_1b<0.01) T_1b=0.01;
  }
  else {
    T_1 = t1;
    T_1b = t1b;
  }

  /*if (inferart) {
    RR = exp( paramcpy(R_index()) );
    if (RR<1) RR=1;
    }*/

  if (inferwm) {
    fwm=paramcpy(wm_index());
    //fwm=20;
    deltwm=paramcpy(wm_index()+1);

    if (infertau) {
      tauwmset = paramcpy(wm_index()+2);
    }
    else tauwmset = seqtau;

    if (infert1) {
      T_1wm = paramcpy(wm_index()+3);
      if (T_1<0.01) T_1=0.01;
    }
    else T_1wm = t1wm;

    if (usepve) {
      pv_gm = paramcpy(pv_index());
      pv_wm = paramcpy(pv_index()+1);
    }
    else {
      pv_gm=1;pv_wm=1;
    }
  }
  else {
    fwm=0;
    deltwm=0;
    T_1wm=t1wm;

    pv_gm=1;
    pv_wm=1;
  }

  p=paramcpy(disp_index());
  s = exp( params(disp_index()+1) );

    float lambdagm = 0.98;
    float lambdawm = 0.82;

    float T_1app = 1/( 1/T_1 + 0.01/lambdagm );
    float T_1appwm = 1/( 1/T_1wm + 0.01/lambdawm );
    //float R = 1/T_1app - 1/T_1b;
    //float Rwm = 1/T_1appwm - 1/T_1b;

    float tau; //bolus length as seen by kintic curve
    float taub; //bolus length of blood as seen in signal
    float tauwm;
    

    //float F=0;
    //float Fwm=0;
    //float Fblood=0;
    float kctissue;
    float kcblood;
    float kcwm;

    //useful time indepedenent terms
    //float gammaps = gamma(1+p*s);
    float k=1+p*s;
    float A = T_1app - T_1b;
    float B = A + s*T_1app*T_1b;
    float C = pow(s-1/T_1app+1/T_1b,p*s);
    float Awm = T_1appwm - T_1b;
    float Bwm = Awm + s*T_1appwm*T_1b;
    float Cwm = pow(s-1/T_1appwm+1/T_1b,p*s);


    // loop over tis
    float ti;
    result.ReSize(tis.Nrows()*repeats);

    for(int it=1; it<=tis.Nrows(); it++)
      {
	ti = tis(it);
	//Fblood = 2*fblood; //scale aBV correctly
	//	F = 2*ftiss;
	//Fwm = 2*fwm;



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

	if(tauwmset < ti -  pretisat)
	  {tauwm = tauwmset; }
	else
	  {tauwm = ti -  pretisat; }


	// (1) tissue contribution 
	    if(ti < delttiss)
	      { kctissue = 0;}
	    else if(ti >= delttiss && ti <= (delttiss + tau))
	      {
		kctissue = 2*ftiss* 1/A * exp( -(T_1app*delttiss + (T_1app+T_1b)*ti)/(T_1app*T_1b) )*T_1app*T_1b*pow(B,-k)*
		  (  exp(delttiss/T_1app + ti/T_1b) * pow(s*T_1app*T_1b,k) * ( 1 - igamc(k,B/(T_1app*T_1b)*(ti-delttiss)) ) +
		     exp(delttiss/T_1b + ti/T_1app) * pow(B,k) * ( -1 + igamc(k,s*(ti-delttiss)) )  );

	      }
	    else //(ti > delttiss + tau)
	      {
		kctissue = 2*ftiss* 1/(A*B) *
		  (  exp(-A/(T_1app*T_1b)*(delttiss+tau) - ti/T_1app)*T_1app*T_1b/C*
		     (  pow(s,k)*T_1app*T_1b* ( -1 + exp( (-1/T_1app+1/T_1b)*tau )*( 1 - igamc(k,B/(T_1app*T_1b)*(ti-delttiss)) ) +
						igamc(k,B/(T_1app*T_1b)*(ti-delttiss-tau)) ) - 
			exp( -A/(T_1app*T_1b)*(ti-delttiss-tau) )*C*B * ( igamc(k,s*(ti-delttiss-tau)) - igamc(k,s*(ti-delttiss)) ) )  );
	      }
	
	    // (2) arterial contribution
	    if(ti < deltblood)
	      { 
		kcblood = 0;
	      }
	    else if(ti >= deltblood && ti <= (deltblood + taub))
	      { 
		kcblood = 2 * fblood * exp(-ti/T_1b) * ( 1 - igamc(k,s*(ti-deltblood))  ); 
	      }
	    else //(ti > deltblood + taub)
	      {
		kcblood = 2 * fblood * exp(-ti/T_1b) * ( igamc(k,s*(ti-deltblood-taub)) - igamc(k,s*(ti-deltblood)) ) ; 
									   
	      }


	    // (3) WM contribution
	    if(ti < deltwm)
	      { kcwm = 0;}
	    else if(ti >= deltwm && ti <= (deltwm + tauwm))
	      {
		kcwm = 2*fwm* 1/Awm * exp( -(T_1appwm*deltwm + (T_1appwm+T_1b)*ti)/(T_1appwm*T_1b) )*T_1appwm*T_1b*pow(Bwm,-k)*
		  (  exp(deltwm/T_1appwm + ti/T_1b) * pow(s*T_1appwm*T_1b,k) * ( 1 - igamc(k,Bwm/(T_1appwm*T_1b)*(ti-deltwm)) ) +
		     exp(deltwm/T_1b + ti/T_1app) * pow(Bwm,k) * ( -1 + igamc(k,s*(ti-deltwm)) )  );

	      }
	    else //(ti > delttiss + tau)
	      {
		kcwm = 2*fwm* 1/(Awm*Bwm) *
		  (  exp(-Awm/(T_1appwm*T_1b)*(deltwm+tauwm) - ti/T_1appwm)*T_1appwm*T_1b/Cwm*
		     (  pow(s,k)*T_1appwm*T_1b* ( -1 + exp( (-1/T_1appwm+1/T_1b)*tauwm )*( 1 - igamc(k,Bwm/(T_1appwm*T_1b)*(ti-deltwm)) ) +
						igamc(k,Bwm/(T_1appwm*T_1b)*(ti-deltwm-tauwm)) ) - 
			exp( -Awm/(T_1appwm*T_1b)*(ti-deltwm-tauwm) )*Cwm*Bwm * ( igamc(k,s*(ti-deltwm-tauwm)) - igamc(k,s*(ti-deltwm)) ) )  );
	      }



	    if (isnan(kctissue)) { kctissue=0; LOG << "Warning NaN in tissue curve at TI:" << ti << " with f:" << ftiss << " delt:" << delttiss << " tau:" << tau << " T1:" << T_1 << " T1b:" << T_1b << " p:" << p << " s:" << s << endl; }
	    if (isnan(kcwm)) { kcwm=0; LOG << "Warning NaN in WM curve at TI:" << ti << " with f:" << fwm << " delt:" << deltwm << " tau:" << tauwm << " T1wm:" << T_1wm << " T1b:" << T_1b << " p:" << p << " s:" << s << endl; }
	    //}

	/* output */
	// loop over the repeats
	for (int rpt=1; rpt<=repeats; rpt++)
	  {
	    result( (it-1)*repeats+rpt ) = pv_gm*kctissue + kcblood + pv_wm*kcwm;
	  }

 
      }
    //cout << result.t();
    

  return;
}


DevelFwdModel::DevelFwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
      repeats = convertTo<int>(args.Read("repeats")); // number of repeats in data
      t1 = convertTo<double>(args.ReadWithDefault("t1","1.3"));
      t1b = convertTo<double>(args.ReadWithDefault("t1b","1.5"));
      t1wm = convertTo<double>(args.ReadWithDefault("t1wm","1.1"));
      lambda =convertTo<double>(args.ReadWithDefault("lambda","0.9")); //NOTE that this parameter is not used!!

      pretisat = convertTo<double>(args.ReadWithDefault("pretisat","0")); // deal with saturation of the bolus a fixed time pre TI measurement
      grase = args.ReadBool("grase"); // DEPRECEATED data has come from the GRASE-ASL sequence - therefore apply pretisat of 0.1s
      if (grase) pretisat=0.1;

      infertau = args.ReadBool("infertau"); // infer on bolus length?
      infert1 = args.ReadBool("infert1"); //infer on T1 values?
      inferart = args.ReadBool("inferart"); //infer on arterial compartment?
      inferwm = args.ReadBool("inferwm");
      //inferinveff = args.ReadBool("inferinveff"); //infer on a linear decrease in inversion efficiency?
      //infertrailing = args.ReadBool("infertrailing"); //infers a trailing edge bolus slope using new model
      seqtau = convertTo<double>(args.ReadWithDefault("tau","1000")); //bolus length as set by sequence (default of 1000 is effectively infinite
      bool ardoff = false;
      ardoff = args.ReadBool("ardoff");
      bool tauboff = false;
      tauboff = args.ReadBool("tauboff"); //forces the inference of arterial bolus off

      usepve = args.ReadBool("usepve");

      // combination options
      infertaub = false;
      if (inferart && infertau && !tauboff) infertaub = true;

      
      //special - turn off tissue cpt
      infertiss=true;
      bool tissoff = args.ReadBool("tissoff");
      if (tissoff) infertiss = false;

      
      // deal with ARD selection
      doard=false;
      tissard=false;artard=true;wmard=true; //default ARD flags
      //if (inferart==true && ardoff==false) { doard=true;}
      //if (inferwm==true && ardoff==false) {doard=true; }
      //special, individual ARD switches
      bool tissardon = args.ReadBool("tissardon");
      if (tissardon) tissard=true;
      bool artardoff = args.ReadBool("artardoff");
      if (artardoff) artard=false;
      bool wmardoff = args.ReadBool("wmardoff");
      if (wmardoff) wmard=false;

      // ** ardoff overrides all other ARD options
      if ( (tissard || artard || wmard) && !ardoff) doard = true;

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
      /*if (tis.Nrows()==1) {
	//only one TI therefore only infer on CBF and ignore other inference options
	LOG << "--Single inversion time mode--" << endl;
	LOG << "Only a sinlge inversion time has been supplied," << endl;
	LOG << "Therefore only tissue perfusion will be inferred." << endl;
	LOG << "-----" << endl;
	singleti = true;
	// force other inference options to be false
	infertau = false; infert1 = false; inferart = false; //inferinveff = false;
	}*/
	
      // add information about the parameters to the log
      LOG << "Inference using development model" << endl;
      if (pretisat>0) LOG << "Saturation of" << pretisat << "s before TI has been specified" << endl;
      if (grase) LOG << "Using pre TI saturation of 0.1 for GRASE-ASL sequence" << endl;
      LOG << "    Data parameters: #repeats = " << repeats << ", t1 = " << t1 << ", t1b = " << t1b;
      LOG << ", bolus length (tau) = " << seqtau << endl ;
      if (infertau) {
	LOG << "Infering on bolus length " << endl; }
      if (doard) {
	LOG << "ARD subsystem is enabled" << endl; }
      if (infertiss) {
	LOG << "Infertting on tissue component " << endl; }
      if (doard && tissard) {
	LOG << "ARD has been set on the tissue component " << endl; }
      if (inferart) {
	LOG << "Infering on artertial compartment " << endl; }
      if (doard && artard) {
	LOG << "ARD has been set on arterial compartment " << endl; }
      if (inferwm) {
	LOG << "Inferring on white matter component" << endl; 
	if (doard && wmard) { LOG << "ARD has been set on wm component" << endl;}
      }
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

void DevelFwdModel::ModelUsage()
{ 
  cout << "\nUsage info for --model=grase:\n"
       << "Required parameters:\n"
       << "--repeats=<no. repeats in data>\n"
       << "--ti1=<first_inversion_time_in_seconds>\n"
       << "--ti2=<second_inversion_time>, etc...\n"
       << "Optional arguments:\n"
       << "--grase *DEPRECEATAED* (data collected using GRASE-ASL: same as --pretissat=0.1)"
       << "--pretisat=<presat_time> (Define that blood is saturated a specific time before TI image acquired)"
       << "--tau=<temporal_bolus_length> (default 10s if --infertau not set)\n"
       << "--t1=<T1_of_tissue> (default 1.3)\n"
       << "--t1b=<T1_of_blood> (default 1.5)\n"
       << "--infertau (to infer on bolus length)\n"
       << "--inferart (to infer on arterial compartment)\n"
       << "--infert1 (to infer on T1 values)\n"
    ;
}

void DevelFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void DevelFwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  if (infertiss) {
  names.push_back("ftiss");
  //if (!singleti) 
    names.push_back("delttiss");
  }
  if (infertau && infertiss)
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
  /*if (inferart) {
    names.push_back("R");
    }*/

  if (inferwm) {
    names.push_back("fwm");
    names.push_back("deltwm");

    if (infertau) names.push_back("tauwm");
    if (infert1) names.push_back("T_1wm");

    if (usepve) {
      names.push_back("p_gm");
      names.push_back("p_wm");
    }
  }
  names.push_back("p");
  names.push_back("s_log");
}

void DevelFwdModel::SetupARD( const MVNDist& theta, MVNDist& thetaPrior, double& Fard)
{
  Tracer_Plus tr("DevelFwdModel::SetupARD");

  if (doard)
    {
      //sort out ARD indices
      if (tissard) ard_index.push_back(tiss_index());
      if (artard) ard_index.push_back(art_index());
      if (wmard) ard_index.push_back(wm_index());

      Fard = 0;

      int ardindex;
      for (unsigned int i=0; i<ard_index.size(); i++) {
	//iterate over all ARD parameters
	ardindex = ard_index[i];

	SymmetricMatrix PriorPrec;
	PriorPrec = thetaPrior.GetPrecisions();
	
	PriorPrec(ardindex,ardindex) = 1e-12; //set prior to be initally non-informative
	
	thetaPrior.SetPrecisions(PriorPrec);
	
	thetaPrior.means(ardindex)=0;
	
	//set the Free energy contribution from ARD term
	SymmetricMatrix PostCov = theta.GetCovariance();
	double b = 2/(theta.means(ardindex)*theta.means(ardindex) + PostCov(ardindex,ardindex));
	Fard += -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
      }
  }
  return;
}

void DevelFwdModel::UpdateARD(
				const MVNDist& theta,
				MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("DevelFwdModel::UpdateARD");
  
  if (doard)
    Fard=0;
    {
      int ardindex;
      for (unsigned int i=0; i<ard_index.size(); i++) {
	//iterate over all ARD parameters
	ardindex = ard_index[i];

  
      SymmetricMatrix PriorCov;
      SymmetricMatrix PostCov;
      PriorCov = thetaPrior.GetCovariance();
      PostCov = theta.GetCovariance();

      PriorCov(ardindex,ardindex) = theta.means(ardindex)*theta.means(ardindex) + PostCov(ardindex,ardindex);

      
      thetaPrior.SetCovariance(PriorCov);

      //Calculate the extra terms for the free energy
      double b = 2/(theta.means(ardindex)*theta.means(ardindex) + PostCov(ardindex,ardindex));
      Fard += -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
    }
  }

  return;

  }

float DevelFwdModel::icgf(float a, float x) const {
  Tracer_Plus("DevelFwdModel::icgf");

  //incomplete gamma function with a=k, based on the incomplete gamma integral

  return gamma(a)*igamc(a,x);
}
