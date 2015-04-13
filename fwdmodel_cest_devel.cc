/*  fwdmodel_cest_devel.cc - Developement CEST APT model

    Michael Chappell, IBME PUMMA & FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_cest_devel.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
#include "miscmaths/miscprob.h"
using namespace NEWIMAGE;
#include "easylog.h"

string CESTDevelFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_cest_devel.cc,v 1.4 2010/10/01 13:46:16 chappell Exp $";
}

void CESTDevelFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("CESTDevelFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors

     // M0
     prior.means(1) = 0.0;
     prior.means(2) = 0.0; //0.001
     prior.means(3) = 0.0;  //0.1;

     precisions(1,1) = 1e-12;
     precisions(2,2) = 1e6; //1e7
     precisions(3,3) = 1e99; //2500; !OFF
     if (mton) precisions(3,3) = 2500;

     // exchnage consts (these are log_e)
     prior.means(4) = 3; //1.3;
     precisions(4,4) = 1;
     prior.means(5) = 3.7;//1.6;
     precisions(5,5) = 10;

     // frequency offsets (ppm)
     prior.means(6) = 0;
     precisions(6,6) = 100;
     
     prior.means(7) = ppm_apt_set; //3.5;
     precisions(7,7) = 100;
     prior.means(8) = -2.34;
     precisions(8,8) = 1e99; // !OFF
     if (mton) precisions(8,8) = 100;
     

     // B1 offset (fractional)
     prior.means(9) = 0;
     precisions(9,9) = 1e14;


     int idx = 10; //should be next entry in priors
     if (pvcorr) {
       // M0
       prior.means(idx) = 0.0;
       prior.means(idx+1) = 0.0;
       prior.means(idx+2) = 0.15;
       
       precisions(idx,idx) = 1e-12;
       precisions(idx+1,idx+1) = 100;
       precisions(idx+2,idx+2) = 2500; 
       
       // exchnage consts (these are log_e)
       prior.means(idx+3) = 3;
       precisions(idx+3,idx+3) = 1;
       prior.means(idx+4) = 3.7;
       precisions(idx+4,idx+4) = 1;

       prior.means(idx+5) = 0.0;
       precisions(idx+5,idx+5) = 1e-12;

       prior.means(idx+6) = 1.0;
       precisions(idx+6,idx+6) = 1e12;
       prior.means(idx+7) = 0.0;
       precisions(idx+7,idx+7) = 1e12;
       prior.means(idx+8) = 0.0;
       precisions(idx+8,idx+8) = 1e12;
       idx += 9;
     }

     //T1 and T2 values
     int t12idx=idx;
     if (t12soft) {
     for (int i=0; i<npool; i++) {
       prior.means(idx+i) = T12master(1,i+1);
       prior.means(idx+npool+i) = T12master(2,i+1);
     }
     
       //T1
       precisions(idx,idx) = 44.4;
       precisions(idx+1,idx+1) = 100;
       precisions(idx+2,idx+2) = 4;
       //T2
       precisions(idx+3,idx+3) = 1e4;
       precisions(idx+4,idx+4) = 4e4;
       precisions(idx+5,idx+5) = 1e12;
       idx += 6;
     }
     //else { //no T12 inference
     //  for (int i=0; i<npool; i++) {
     //	 precisions(idx+i,idx+i) = 1e99;
     //  }
     //}

     int wassridx=idx;
     if (wassron & !wassronly) {
       //M0 for WASSR image
       prior.means(idx) = 0;
       precisions(idx,idx) = 1e-12;
       idx +=1;
      }

     // --- Special cases ---
     if (wassronly) {
       //special case doing wassr only
       precisions = IdentityMatrix(NumParams()) * 1e99; //set the precisions to all known
       precisions(1,1) = 1e-12; //release the precision on M0 water
       precisions(6,6)=100; // release the precision on the w offset

       if (t12soft) {
	 precisions(t12idx,t12idx) = 44.4; //release the precision on T1 (water)
	 precisions(t12idx+3,t12idx+3) = 1e4; // precision on T2 (water)
       }
     } 

     if (basic) {
       //special case basic water only fit
        precisions = IdentityMatrix(NumParams()) * 1e99; //set the precisions to all known
	precisions(1,1) = 1e-12; //release the precision on M0 water
	precisions(6,6)=100; // release the precision on the w offset
	precisions(9,9)=1e14; //release the precisions for B1

	prior.means(3)=0.0; //turn off MT

	if (t12soft) {
	  precisions(t12idx,t12idx) = 44.4; //release the precision on T1 (water)
	  precisions(t12idx+3,t12idx+3) = 1e4; // precision on T2 (water)
       }
     }

    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital posterior
      posterior.means(1) = 100;
      precisions(1,1) = 1e-4;
    if (t12soft | wassronly) {
    }
    else {
      posterior.means(2) = 0.001;
      precisions(2,2) = 1e7;
    }

    if (wassron & !wassronly) {
       //M0 for WASSR image
       prior.means(wassridx) = 100;
       precisions(wassridx,wassridx) = 1e-4;
      }


    posterior.SetPrecisions(precisions);
    
}    
    
    

void CESTDevelFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("CESTDevelFwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
   for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
      }
  

  //ColumnVector wi_fixed(npool);
  //wi_fixed << 42.58e6*3  << 42.58e6*3*3.5e-6;
  //wi_fixed << wlam << wlam + ppm_apt*wlam/1e6;

   //model matrices
   ColumnVector M0(npool);
   Matrix kij(npool,npool);
   Matrix T12(2,npool);

   // extract values from params
   // M0 comes first
   int place=1;
   M0(1) = paramcpy(place); // this is the 'master' M0 value of water
   if (M0(1)<1e-4) M0(1)=1e-4; //M0 of water cannot disapear all together
   place++;
   // this code here if parameters contains actual M0 values
   //M0.Rows(2,npool) = paramcpy.Rows(place,place+npool-1-1);
   // place +=npool-1;

   // values in the parameters are ratios of M0_water
   float M0ratio;
   for (int j=2; j<=npool; j++) {
     M0ratio = paramcpy(place);
     if (M0ratio > 0.1) { //dont expect large ratios
       M0ratio = 0.1;
     }
     M0(j) = M0ratio*M0(1);

     place++;
     }

   // now exchange - we assume that only significant exchnage occurs with water
   kij=0.0; //float ktemp;
   for (int j=2; j<=npool; j++) {
     //ktemp = exp(params(place)); //this is the 'fundamental rate const - non-linear transformation
     //kij(j,1) = ktemp*M0(1);
     //kij(1,j) = ktemp*M0(j);
     kij(j,1) = exp(params(place));  //non-linear transformation
     kij(1,j) =kij(j,1)*M0(j)/M0(1);
     
     place++;
     
   }

   // frequency offset next
   float ppm_off = params(place); //(params.Rows(place,place)).AsScalar(); // frequncy offset due to field
    place++;
    float ppm_apt = params(place);
    place++;
    float ppm_mt = params(place);
    place++;

    // now B1 offset
    float B1off = params(place)*1e6; //scale the B1_off parameter to achieve proper updating of this parameter
                                     // (otherwise it is sufficently small that numerical diff is ineffective)
    place++;

    // PV correction section (untested)
	ColumnVector M0_WM(npool);
	Matrix kij_WM(npool,npool);
	Matrix T12_WM(2,npool);
	Matrix T12_CSF(2,1);
	Matrix M0_CSF(1,1);
	float PV_GM;
	float PV_WM;
	float PV_CSF;
    if (pvcorr) 
      {

	// now parameters for WM
	M0_WM(1) = paramcpy(place); // this is the M0 value of water
	if (M0_WM(1)<1e-4) M0_WM(1)=1e-4; //M0 of water cannot disapear all together
	place++;
	// this code here if parameters contains actual M0 values
	//M0.Rows(2,npool) = paramcpy.Rows(place,place+npool-1-1);
	// place +=npool-1;
	
	// values in the parameters are ratios of M0_water
	float M0WMratio;
	for (int j=2; j<=npool; j++) {
	  M0WMratio = paramcpy(place);
	  if (M0WMratio > 0.1) { //dont expect large ratios
	    M0WMratio = 0.1;
	  }
	  M0_WM(j) = M0WMratio*M0_WM(1);
	  
	  place++;
	}
	
	// now exchange - we assume that only significant exchnage occurs with water
	kij_WM=0.0; //float ktemp;
	for (int j=2; j<=npool; j++) {
	  //ktemp = exp(params(place)); //this is the 'fundamental rate const - non-linear transformation
	  //kij(j,1) = ktemp*M0(1);
	  //kij(1,j) = ktemp*M0(j);
	  kij_WM(j,1) = exp(params(place));  //non-linear transformation
	  kij_WM(1,j) =kij_WM(j,1)*M0_WM(j)/M0_WM(1);
	  
	  place++;
	}

	  // CSF
	  M0_CSF = paramcpy(place);
	  place++;

	  //partial volume estimates
	  PV_GM = paramcpy(place);
	  place++;
	  PV_WM = paramcpy(place);
	  place++;
	  PV_CSF = paramcpy(place);
	  place++;
	  
	  //T12 fixed
	  T12_WM = T12WMmaster;
	  T12_CSF = T12CSFmaster;
      }

    // T12 parameter values
    if (t12soft) {
      // T12 values
      for (int i=1; i<=npool; i++) {
	T12(1,i) = paramcpy(place);
	if (T12(1,i)<1e-12) T12(1,i)=1e-12; // 0 is no good for a T1 value
	if (T12(1,i)>10) T12(1,i)=10; // Prevent convergence issues causing T1 to blow up
	place++;
      }
      for (int i=1; i<=npool; i++) {
	T12(2,i) = paramcpy(place);
	if (T12(2,i)<1e-12) T12(2,i)=1e-12; //0 is no good for a T2 value

	if (T12(2,i)>1) T12(2,i)=1; // Prevent convergence issues causing T2 to blow up
	place++;
      }
       }
         else {T12 = T12master;}

    //WASSR image M0
    float M0_WASSR=0.0;
    if (wassron & !wassronly) {
      M0_WASSR = paramcpy(place);
      place++;
    }
   else {
     M0_WASSR = (M0.Row(1)).AsScalar(); //if we are doing wassr only then use the main M0 value
     if (wassron) place++;
   }
      

   //cout << "Parameters set up" << endl;
   //cout << "M0: " << M0.t() << endl << "wi: " << wi.t() << endl << "kij: " << kij << endl;

   //MODEL - CEST
    ColumnVector cest_result;

    //no saturation image first
    ColumnVector nosat(1);
    //cest_result.ReSize(1);
    if (pvcorr) {
      nosat = PV_GM*M0.Row(1) + PV_WM*M0_WM.Row(1) + PV_CSF*M0_CSF.Row(1);}
    else {
      nosat = M0.Row(1);}
    

   if (!wassronly) {

     //iterate over the data sets
     for (unsigned int n=0; n<t.size(); n++)
       {
	 // deal with frequencies
	 ColumnVector wi(npool);
	 float wlocal = wlam[n]*ppm_off/1e6; //local water frequency
	 //cout << wlocal << endl;
	 wi << wlocal << wlam[n]*ppm_apt/1e6 + wlocal*(1+ppm_apt/1e6) << wlam[n]*ppm_mt/1e6 + wlocal*(1+ppm_mt/1e6); // species b is at ppm*wlam, but alos include offset of main field
	 //cout << wi << endl;

	 //deal with B1
	 if (B1off<-0.5) B1off=-0.5; // B1 cannot go too small
	 if (B1off>10) B1off=10; //unlikely to get this big (hardlimit in case of convergence problems)
	 float B1 = B1set[n] * (1+B1off);
	 float w1 = 42.58e6*B1*2*M_PI; // in radians!

	 ColumnVector thisresult;
	 if (pvcorr) {
	   ColumnVector GM_result;
	   GM_result = Mz_spectrum(wvec[n],w1,t[n],M0,wi,kij,T12);
	   ColumnVector WM_result;
	   WM_result = Mz_spectrum(wvec[n],w1,t[n],M0_WM,wi,kij_WM,T12_WM);
	   ColumnVector CSF_result;
	   CSF_result = Mz_spectrum(wvec[n],w1,t[n],M0_CSF,wi.Row(1),kij.SubMatrix(1,1,1,1),T12_CSF);

	   thisresult = PV_GM*GM_result + PV_WM*WM_result + PV_CSF*CSF_result;
	 }
	 else {
	   thisresult = Mz_spectrum(wvec[n],w1,t[n],M0,wi,kij,T12);
	 }

	 cest_result &= thisresult;
	 /*
	 if (n>0) {
	   cest_result = cest_result & thisresult;
	 }
	 else { cest_result = thisresult; }
	 */
       }

     cest_result &= nosat;
   }

   //cout << cest_result.t() << endl;

 //MODEL - WASSR
   ColumnVector wassr_result;
   
   //no saturation image first
   ColumnVector wassr_nosat(1);
   //wassr_result.ReSize(1);
   wassr_nosat = M0_WASSR;
   

   if (wassron) {
     //wassr_result = Mz_spectrum(W_wvec,W_w1,W_t,M0,wi,kij,T12); //NB this line uses the main CEST wi not specific to the WASSR acquisition centre frequency
     Matrix W_wlocal(1,1);
     ColumnVector M0W(1);
     M0W = M0_WASSR;
     W_wlocal = W_wlam*ppm_off/1e6;
     wassr_result &= Mz_spectrum(W_wvec,W_w1,W_t,M0W,W_wlocal,kij.SubMatrix(1,1,1,1),T12.Column(1)); //use simple 1-pool model for the wassr spectrum
   }

   wassr_result &= wassr_nosat;
     
   if (wassronly) { result = wassr_result; }
   else { 
     if(wassron) { result = cest_result & wassr_result;}
     else        { result = cest_result; }
   }


  return;
}


CESTDevelFwdModel::CESTDevelFwdModel(ArgsType& args)
{
  Tracer_Plus tr("CESTDevelFwdModel");

    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    
    wassron=false; t12soft=false; pvcorr=false;basic=false;mton=false;

    // define input parameters that we need here 
    vector<float> wrange;
    vector<float> wstart;
    vector<float> wstep;
    vector<int> nfreq;

    //WASSR input parameters
    float W_B1=0.0; 
    float W_wrange;
    int W_nfreq;

    vector<float> wcentre;
    float W_wcentre=-1.0;

    bool revfrq;

    if (scanParams == "cmdline")
    {
      //primary dataset
      string tsat_text = args.ReadWithDefault("tsat","none!");
      if (tsat_text != "none!") {
	B1set.push_back( convertTo<double>(args.ReadWithDefault("B1","0.0")) );
	t.push_back( convertTo<double>(tsat_text) );
	wrange.push_back( convertTo<double>(args.ReadWithDefault("wrange","0.0")) );
	wstart.push_back( convertTo<double>(args.ReadWithDefault("wstart","0.0")) );
	wstep.push_back( convertTo<double>(args.ReadWithDefault("wstep","0.0")) );
	nfreq.push_back( convertTo<int>(args.ReadWithDefault("nfreq","0")) );
	wcentre.push_back( convertTo<double>(args.ReadWithDefault("wcentre","-1.0")) );
      }

      //additional datasets
      int N=1;
      while (true) 
	{
	  tsat_text = args.ReadWithDefault("tsat_"+stringify(N),"stop!");
	  //cout << tsat_text << endl;
	  if (tsat_text == "stop!") break; //there are no more datasets defined

	  //append new dataset parameters
	  t.push_back( convertTo<double>(tsat_text) );
	  B1set.push_back( convertTo<double>(args.ReadWithDefault("B1_"+stringify(N),"0.0")) );
	  wrange.push_back( convertTo<double>(args.ReadWithDefault("wrange_"+stringify(N),"0.0")) );
	  wstart.push_back( convertTo<double>(args.ReadWithDefault("wstart_"+stringify(N),"0.0")) );
	  wstep.push_back( convertTo<double>(args.ReadWithDefault("wstep_"+stringify(N),"0.0")) );
	  nfreq.push_back( convertTo<int>(args.ReadWithDefault("nfreq_"+stringify(N),"0")) );
	  wcentre.push_back( convertTo<double>(args.ReadWithDefault("wcentre_"+stringify(N),"-1.0")) );
	  N++;
	}

	  
      W_B1 = convertTo<double>(args.ReadWithDefault("W_B1","0.0"));
      W_t = convertTo<double>(args.ReadWithDefault("W_tsat","0.0"));
      W_wrange = convertTo<double>(args.ReadWithDefault("W_wrange","0.0"));
      W_nfreq = convertTo<int>(args.ReadWithDefault("W_nfreq","0"));
      W_wcentre = convertTo<double>(args.ReadWithDefault("W_wcentre","-1.0"));

      wassronly = args.ReadBool("wassronly");
      basic  = args.ReadBool("basic");

      pvcorr = args.ReadBool("pvcorr");

      t12soft = args.ReadBool("t12prior");

      mton = args.ReadBool("infermt");

      ppm_apt_set = convertTo<double>(args.ReadWithDefault("ppm_apt","3.5"));

      revfrq=args.ReadBool("reverse");
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    

    cout << "No. CEST datasets: " << t.size() << endl;

    //for (unsigned int n=0; n<t.size(); n++)
    //  { cout << t[n] << endl;
    //	cout << wcentre[n] << endl;}

    //determine if we are doing WASSR
    if (W_B1>1e-12) {
      wassron=true;
      //cout << "WASSR is on" << endl;
    }

    //initialization
    npool = 3;

    // some fixed things
    T12master.ReSize(2,npool);
    T12master << 1.3 << 0.77  << 1  
	      << 0.05 << 0.01 << 1e-5;  //T2

    T12WMmaster.ReSize(2,npool);
    T12WMmaster << 1 << 0.77  << 1  
	      << 0.02 << 0.01 << 1e-5;  //T2
    T12CSFmaster.ReSize(2,1);
    T12CSFmaster << 3.5 << .75;
    
    //ppm_apt = -3.5;
    //ppm_mt = -2.41;

    //deal with water centre frequency *in radians!*
    float wdefault = 42.58e6*3*2*M_PI; // this is the default values (3T)
    if (W_wcentre>0)   W_wlam = W_wcentre*1e6*2*M_PI; //NB input is in MHz
    else                W_wlam = wdefault;

    wlam.resize(t.size());
    for (unsigned int n=0; n<t.size(); n++) {
      //iterate over all the datasets

      //water centre frequency
      if (wcentre[n]>0)     wlam[n] = wcentre[n]*1e6*2*M_PI; //NB input is in MHz
      else                  wlam[n] = wdefault;

      // CEST frequency vector
      //int nfreq=32;
      ColumnVector ppmvec(nfreq[n]);
      float ppminc;
      if (wrange[n] > 0.0) 
	{ //deafult case is that wrange is specified
	  ppmvec(1) = -wrange[n];
	  ppminc = (2*wrange[n])/(nfreq[n]-1);
	}
      else
	{ //if wrange is zero then we use wstart and wstep parameters
	  ppmvec(1) = wstart[n];
	  ppminc = wstep[n];
	}
      //cout << ppminc << endl;
      for (int i=2; i<=nfreq[n]; i++) {
	ppmvec(i) = ppmvec(i-1) + ppminc;
      }
   
      //wvec.ReSize(nfreq);
      ColumnVector freqvec(nfreq[n]);
      freqvec = ppmvec*wlam[n]/1e6;

      if (revfrq) freqvec = -freqvec;
      //cout << freqvec.t() << endl;
      wvec.push_back( freqvec );
      //cout << wvec[n] << endl;
      //cout << wlam[n] << endl;
    }


    if (wassron) {
      //float W_B1 = 0.4e-6;
      W_w1 = 42.58e6*W_B1*2*M_PI; // in radians!

      // WASSR frequency vector
      ColumnVector ppmvec(W_nfreq);
      ppmvec(1) = -W_wrange;
      float ppminc = (2*W_wrange)/(W_nfreq-1);
      //cout << ppminc << endl;
      for (int i=2; i<=W_nfreq; i++) {
	ppmvec(i) = ppmvec(i-1) + ppminc;
      }
      W_wvec.ReSize(W_nfreq);
      W_wvec = ppmvec*W_wlam/1e6;
    }
}

void CESTDevelFwdModel::ModelUsage()
{ 
  cout << "\nUsage info for --model=CEST:\n"
       << "Undefined\n"
    ;
}

void CESTDevelFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void CESTDevelFwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  names.push_back("M0a");
  names.push_back("M0b_r");
  names.push_back("M0c_r");
  names.push_back("kba");
  names.push_back("kca");
  names.push_back("ppm_off");
  names.push_back("ppm_apt");
  names.push_back("ppm_mt");
  names.push_back("B1_off");
  if (pvcorr) {
    names.push_back("WM_M0a");
    names.push_back("WM_M0b_r");
    names.push_back("WM_M0c_r");
    names.push_back("WM_kba");
    names.push_back("WM_kca");
    names.push_back("CSF_M0a");
  }
  if (t12soft) {
    names.push_back("T1a");
  names.push_back("T1b");
  names.push_back("T1c");
  names.push_back("T2a");
  names.push_back("T2b");
  names.push_back("T2c");
  }
  if (wassron) {
    names.push_back("M0_WASSR");
  }
}

void CESTDevelFwdModel::SetupARD( const MVNDist& theta, MVNDist& thetaPrior, double& Fard)
{
  Tracer_Plus tr("CESTDevelFwdModel::SetupARD");

  if (doard)
    {
      //sort out ARD indices

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

void CESTDevelFwdModel::UpdateARD(
				const MVNDist& theta,
				MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("CESTDevelFwdModel::UpdateARD");
  
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

/*

ReturnMatrix CESTDevelFwdModel::expm(Matrix inmatrix) const
{
  // Do matrix exponential
  // Algorithm from Higham, SIAM J. Matrix Analysis App. 24(4) 2005, 1179-1193
  Tracer_Plus tr("CESTDevelFwdModel::expm");

  Matrix A = inmatrix;

  //Coeff of degree 13 Pade approximant
  //ColumnVector b(13);
  //  b= PadeCoeffs(13);
  
  Matrix X(A.Nrows(),A.Ncols());
  int metflag=0;
  int mvals[5] = {3, 5, 7, 9, 13};
  //ColumnVector thetam(5);
  float thetam[5] = {1.495585217958292e-002, 2.539398330063230e-001, 9.504178996162932e-001, 2.097847961257068e+000, 5.371920351148152e+000};
  int i=0;

  while (metflag<1)
    {
      if (i<4) {
	//cout << "Try m = " << mvals[i] << endl;
	if (A.Norm1() <= thetam[i]) {
	  X = PadeApproximant(A,mvals[i]);
	  metflag=1;
	}
      }
      else {
	// Using Pade Approximant degree 13 and scaling and squaring
	//cout << "Doing m = 13" << endl;
	int s = ceil(log2(A.Norm1()/thetam[4]));
	//cout << s << endl;
	float half=0.5;
	A *= MISCMATHS::pow(half,s);
	X = PadeApproximant(A,13);
	for (int i=1; i<=s; i++) {
	  X *= X;
	}
	metflag=1;
      }
      i++;
  }

  return X;
}

ReturnMatrix CESTDevelFwdModel::PadeApproximant(Matrix inmatrix, int m) const
{
  Tracer_Plus tr("CESTDevelFwdModel::PadeApproximant");

  //cout << "PadeApproximant" << endl;
  //cout << inmatrix << endl;
  assert(inmatrix.Nrows()==inmatrix.Ncols());
  int n=inmatrix.Nrows();
  IdentityMatrix I(n);
  ColumnVector coeff = PadeCoeffs(m);
  Matrix X(n,n);
  Matrix U(n,n); Matrix V(n,n);
  vector<Matrix> Apowers;

  //cout << coeff << endl;

  switch (m)
    {
    case 3:
    case 5:
    case 7:
    case 9:
      //cout << "case " << m << endl;
      Apowers.push_back(I);
      Apowers.push_back(inmatrix*inmatrix);
      for (int j=3; j<=ceil((m+1/2)); j++) {
	Apowers.push_back(Apowers[j-2]*Apowers[1]);
      }
      U=0.0; V=0.0;
      //cout << "Apowers" <<endl;
      for (int j=m+1; j>=2; j -= 2) {
	U += coeff(j)*Apowers[j/2-1];
      }
      U = inmatrix*U;
      for (int j=m; j>=1; j -= 2) {
	V += coeff(j)*Apowers[(j+1)/2-1];
      }
      //cout << U << endl << V << endl;
      X = (U+V)*(-U+V).i();
      //cout << X<< endl;
      break;
    case 13:
      //cout << "case 13" << endl;
      Matrix A2(n,n); Matrix A4(n,n); Matrix A6(n,n);
      A2 = inmatrix*inmatrix; A4 = A2*A2; A6 = A2*A4;
      U = inmatrix * (A6*(coeff(14)*A6 + coeff(12)*A4 + coeff(10)*A2) + coeff(8)*A6 + coeff(6)*A4 + coeff(4)*A2 + coeff(2)*I );
      //cout << U;
      V = A6*(coeff(13)*A6 + coeff(11)*A4 + coeff(9)*A2) + coeff(7)*A6 + coeff(5)*A4 + coeff(3)*A2 + coeff(1)*I;
      //cout << V;
      X = (U+V)*(-U+V).i();
      //cout << X <<endl;
      break;
    }

  return X;
}

ReturnMatrix CESTDevelFwdModel::PadeCoeffs(int m) const {

  Tracer_Plus tr("CESTDevelFwdModel::PadeCoeffs");
  ColumnVector C;
  C.ReSize(m+1);

  //cout << "PadeCoeffs" << endl;

switch (m)
    {
    case 3:
      C << 120 << 60 << 12 << 1;
      break;
    case 5:
      C << 30240 << 15120 << 3360 << 420 << 30 << 1;
      break;
    case 7:
	C << 17297280 << 8648640 << 1995840 << 277200 << 25200 << 1512 << 56 << 1;
      break;
    case 9:
	C << 1.7643225600e10 << 8.821612800e9 << 2.075673600e9 << 3.02702400e8 << 30270240 << 2162160 << 110880 <<3960 << 90 << 1;
      break;
    case 13:
      C << 6.4764752532480000e16 << 3.2382376266240000e16 << 7.771770303897600e15 << 1.187353796428800e15 <<  1.29060195264000e14 <<   1.0559470521600e13 <<  6.70442572800e11 <<      3.3522128640e10 << 1323241920 << 40840800 << 960960 << 16380 << 182 << 1;
      break;
    }
 return C;
}
*/

ReturnMatrix CESTDevelFwdModel::Mz_spectrum(ColumnVector wvec, float w1, float t, ColumnVector M0, ColumnVector wi, Matrix kij, Matrix T12) const {

  Tracer_Plus tr("CESTDevelFwdModel::Mz_spectrum");
  int nfreq = wvec.Nrows();
  int mpool=M0.Nrows();

  //assmeble model matrices
  ColumnVector k1i(mpool);
   ColumnVector k2i(mpool);
   for (int i=1; i<=mpool; i++) {
     k1i(i) = 1/T12(1,i) + (kij.Row(i)).Sum();
     k2i(i) = 1/T12(2,i) + (kij.Row(i)).Sum();
   }

   //cout << "k matrices generated" << endl;
   //cout << k1i << endl << k2i << endl;

   // first population of A (with w=0)
   Matrix A(mpool*3,mpool*3);
   A = 0.0; int st=0;
   for (int i=1; i<=mpool; i++) {
     Matrix D(3,3);
     D=0.0;
     D(1,1) = -k2i(i); D(2,2) = -k2i(i);
     D(1,2) = -(wi(i)); D(2,1) = wi(i);
     D(2,3) = -w1; D(3,2) = w1;
     D(3,3) = -k1i(i);
     //cout << D << endl;
     st = (i-1)*3;
     A.SubMatrix(st+1,st+3,st+1,st+3) = D;
   }
   //cout << A << endl;

   int st2=0; IdentityMatrix I(3);
   for (int i=1; i<=mpool; i++) {
     for (int j=1; j<=mpool; j++) {
       if (i!=j) {
	 st = (i-1)*3;
	 st2 = (j-1)*3;
	 A.SubMatrix(st+1,st+3,st2+1,st2+3) = I*kij(j,i); //NB 'reversal' of indices is correct here
       }
     }
   }

   //cout << "Inital A matrix populated" << endl;
   //cout << A << endl;

   ColumnVector M0i(mpool*3);
   ColumnVector B(mpool*3);
   B = 0.0; M0i=0.0;
   for (int i=1; i<=mpool; i++) {
     M0i(i*3) = M0(i);
     B(i*3) = M0(i)/T12(1,i);
   }

   Matrix M(mpool*3,nfreq);
   M=0.0;
   Matrix AinvB(mpool,mpool); AinvB=0.0;
   for (int k=1; k<=nfreq; k++) {
     //Calculate new A matrix
     for (int i=1; i<=mpool; i++) {
     st = (i-1)*3;
     A(st+1,st+2) = -(wi(i)-wvec(k));
     A(st+2,st+1) = wi(i)-wvec(k);
   }
     //cout << A << endl;
     AinvB = A.i()*B;
     M.Column(k) = expm(A*t) * (M0i + AinvB) - AinvB;
   }

   ColumnVector result;
   //cout << M.Row(3);
   result = (M.Row(3)).AsColumn();
   return result;
}

