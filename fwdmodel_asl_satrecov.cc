/*  fwdmodel_asl_satrecov.cc - Saturation Recovery curve calibration for ASL

    Michael Chappell, IBME & FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_satrecov.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
#include "miscmaths/miscprob.h"
using namespace NEWIMAGE;
#include "easylog.h"

string SatrecovFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_asl_satrecov.cc,v 1.4 2012/12/19 17:16:43 chappell Exp $";
}

void SatrecovFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("SatrecovFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

     prior.means(1)=0;

     prior.means(2) = t1;
     precisions(2,2) = 10; //1e-12;

     prior.means(3) = 1;
     if (fixA) {
       precisions(3,3) = 1e12;
     }
     else {
       precisions(3,3) = 10;
     }

     if (LFAon) {
       prior.means(4) = 1;
       precisions(4,4) = 100;
     }

    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    
    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital posterior
    posterior.means(1)=10;
    precisions(1,1) = 0.1;
    
    posterior.SetPrecisions(precisions);
    
}    
    
    

void SatrecovFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("SatrecovFwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
   for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
      }
  
   float M0t;
   float T1t;
   float A;
   float FA;
   float lFA;
   float g;

   M0t = paramcpy(1);
   T1t = paramcpy(2);
   A = paramcpy(3);

   if (LFAon) {
     g = paramcpy(4);
   }
   else g=1.0;

   //if (g<0.5) g=0.5;
   //if (g>1.5) g=1.5;

   FA =(g+dg)* FAnom;
   lFA = (g+dg)* LFA;

   float T1tp = T1t;
   float M0tp = M0t;

   if (looklocker) {
     T1tp = 1/( 1/T1t - log(cos(FA))/dti );
     M0tp = M0t*(1 - exp(-dti/T1t) )/(1 - cos(FA)*exp(-dti/T1t));
     // note that we do not have sin(FA) here - we actually estiamte the M0 at the flip angle used for the readout!
   }

    // loop over tis
    //float ti;
   if (LFAon) result.ReSize(tis.Nrows()*(nphases+1)*repeats);
   else result.ReSize(tis.Nrows()*nphases*repeats);

   int nti=tis.Nrows();
   double ti;

    for (int ph=1; ph<=nphases; ph++) {
      for (int it=1; it<=tis.Nrows(); it++) {
	for (int rpt=1; rpt<=repeats; rpt++)
	  {
	    ti = tis(it) + slicedt*coord_z; //account here for an increase in delay between slices
	    result( (ph-1)*(nti*repeats) + (it-1)*repeats+rpt ) = M0tp*(1-A*exp(-ti/T1tp));
	  }
      }
    }
    if (LFAon) {
      int ph=nphases+1;
      T1tp = 1/( 1/T1t - log(cos(lFA))/dti );
      M0tp = M0t*(1 - exp(-dti/T1t) )/(1 - cos(lFA)*exp(-dti/T1t));
      for (int it=1; it<=tis.Nrows(); it++) {
	for (int rpt=1; rpt<=repeats; rpt++)
	  {
	    ti = tis(it) + slicedt*coord_z; //account here for an increase in delay between slices
	    result( (ph-1)*(nti*repeats) + (it-1)*repeats+rpt ) = M0tp*sin(lFA)/sin(FA)*(1-A*exp(-tis(it)/T1tp));
	    //note the sin(LFA)/sin(FA) term since the M0 we estimate is actually MOt*sin(FA)
	  }
      }
    }

    

  return;
}


SatrecovFwdModel::SatrecovFwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
      repeats = convertTo<int>(args.ReadWithDefault("repeats","1")); // number of repeats in data
      t1 = convertTo<double>(args.ReadWithDefault("t1","1.3"));
      nphases = convertTo<int>(args.ReadWithDefault("phases","1"));
      slicedt = convertTo<double>(args.ReadWithDefault("slicedt","0.0")); // increase in TI per slice

      fixA=args.ReadBool("fixa"); //to fix the A parameter where it will be ambiguous

      // with a look locker readout
      FAnom = convertTo<double>(args.ReadWithDefault("FA","0"));
      looklocker=false;
      if (FAnom>0.1) looklocker=true;
      cout << "Looklocker" << looklocker << endl;
      FAnom = FAnom * M_PI/180;
      LFA = convertTo<double>(args.ReadWithDefault("LFA","0"));
      LFA = LFA * M_PI/180;
      LFAon=false;
      if (LFA>0) LFAon=true;
      
      dg=0.023;

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
      dti = tis(2)-tis(1); //assuming even sampling!! - this only applies to LL acquisitions
      
      // need to set the voxel coordinates to a deafult of 0 (for the times we call the model before we start handling data)
      coord_x = 0;
      coord_y = 0;
      coord_z = 0;
  
	  
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    
    
 
}

void SatrecovFwdModel::ModelUsage()
{ 
  cout << "\nUsage info for --model=:\n"
       << "Required parameters:\n"
       << "--ti{n}=<nth_inversion_time_in_seconds>\n"
       << "Optional arguments:\n"
       << "--repeats=<no. repeats in data>  {default:1}\n"
       << "--phases=<no. of phases in data> {default:1}\n"
       << "--t1=<T1_of_tissue_prior_mean> {default 1.3}\n"
       << "--looklocker Data was aquired using Look Locker readout\n"
       << "--fa=<flip_angle_in_degrees> for LL readout\n"
       << "--lfa=<low_flip_angle_in_degrees> extra phase in data with low FA\n"
    ;
}

void SatrecovFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void SatrecovFwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  names.push_back("M0t");
  names.push_back("T1t");
  names.push_back("A");
  if (LFAon) { names.push_back("g"); }
}

