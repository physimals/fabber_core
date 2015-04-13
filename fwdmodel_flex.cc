/*  fwdmodel_flex.cc - FLEX model

    Michael Chappell, IBME PUMMA & FMRIB Image Analysis Group

    Copyright (C) 2011 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_flex.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
#include "miscmaths/miscprob.h"
using namespace NEWIMAGE;
#include "easylog.h"

string FLEXFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_flex.cc,v 1.2 2012/02/09 11:50:09 chappell Exp $";
}

void FLEXFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("FLEXFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors

     int place=1;
     // baseline
     for (int i=1; i<=npoly; i++) {
       prior.means(place) = 0.0;
       precisions(place,place) = 0.01;
       place++;
     }

     // PTR
     for (int s=1; s<=ncomp; s++) {
       prior.means(place) = 0.0;
       place++;
     }
     //NB dont really need to worry about these are they will be ARD

     // deltw (ppm)
     for (int s=1; s<=ncomp; s++) {
       prior.means(place) = compspec(s,1);
       if (inferdw) {
	 precisions(place,place) = 10;
       }
       else { precisions(place,place) = 1e12; }
       place++;
     }

     // kevol (log)
     for (int s=1; s<=ncomp; s++) {
       prior.means(place) = compspec(s,2);
       if (inferk) {
	 precisions(place,place) = 1;
       }
       else { precisions(place,place) = 1e12; }
       place++;
     }

     //phase
     for (int s=1; s<=ncomp; s++) {
       prior.means(place) = 0;
       precisions(place,place) = 1;
       place++;
     }

    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital posterior
      posterior.means(1) = 0;
      precisions(1,1) = 10;
      place=2;

      if (npoly>1) {
	for (int i=1; i<=npoly; i++) {
	  precisions(place,place)=10;
	  place++;
	}
      }

            for (int s=1; s<=ncomp; s++) {
	      posterior.means(place) = 1;
	      precisions(place,place) = 10;
	      place++;
            }

      posterior.SetPrecisions(precisions);
    
}    
    
    

void FLEXFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("FLEXFwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
  ColumnVector paramcpy = params;
   for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
      }

   //model matrices
   ColumnVector baseline(npoly);
   ColumnVector PTR(ncomp);
   ColumnVector deltw(ncomp);
   ColumnVector kevol(ncomp);
   ColumnVector phase(ncomp);

   // extract values from params
   int place=1;
   baseline = params.Rows(place,place+npoly-1);
   place += npoly;

   // PTR
   PTR = params.Rows(place,place+ncomp-1)/1e3;
   place += ncomp;

   // deltw
   deltw = params.Rows(place,place+ncomp-1);
   //ColumnVector pw;
   //pw = params.Rows(place,place+ncomp-1);
   //for (int i=1; i<=ncomp; i++) {
     //if (pw(i)>M_PI/2-1e-12) pw(i) = M_PI/2-1e-12;
     //if (pw(i)<-M_PI/2+1e-12) pw(i) = -M_PI/2+1e-12;
     //deltw(i) = compspec(i,1) + tan(pw(i));
   //  deltw(i) = compspec(i,1) + pw(i);
   //}
   
   //cout << deltw.t() << endl;

   place += ncomp;
   deltw /= 1e6; // starts out in ppm
   deltw *= field; //into Hz
   deltw = o1 - deltw; //offset from o1 frequency
   deltw *= 2*M_PI; //into radians/s

   // kevol = (ksw-1/T*_2s)
   kevol = exp( params.Rows(place,place+ncomp-1) ); //infer log(kevol)
   place += ncomp;

   //phase
   phase = params.Rows(place,place+ncomp-1);
   place += ncomp;


    result.ReSize(ntpts);
    result = 0.0;
    
    //baseline
    ColumnVector tpower(tevol);
    tpower=1;
    for (int i=1; i<=npoly; i++) {
      result += baseline(i)*tpower;
      tpower = SP(tpower,tevol);
    }

    //cout << deltw.t() << endl;

    for (int s=1; s<=ncomp; s++) {
      for (int i=1; i<=ntpts; i++) {
	result.Row(i) += PTR.Row(s)*exp( -kevol.Row(s)*tevol.Row(i) )*cos( (deltw.Row(s)*tevol.Row(i) + phase.Row(s)).AsScalar());
      }
    }

    if (ffton) {
      ColumnVector fftx;
      ColumnVector ffty;
      RealFFT(result,fftx,ffty);
      result = sqrt(SP(fftx,fftx)+SP(ffty,ffty));
    }

    //cout << params.t() << endl;
    //cout << result.t() << endl;

  return;
}


FLEXFwdModel::FLEXFwdModel(ArgsType& args)
{
  Tracer_Plus tr("FLEXFwdModel");

    string scanParams = args.ReadWithDefault("scan-params","cmdline");
      bool ardon;

    // compspec
    // 2 Columns: Freq (ppm), kevol (s^-1)
    // Nrows = number of components

    if (scanParams == "cmdline")
    {
      // read timings from file
      tevol = read_ascii_matrix(args.Read("tevol")); //in s

      //read component specificaiton from file
      compspec = read_ascii_matrix(args.Read("comps"));

      //read field (Hz)
      field = convertTo<double>(args.Read("field"));

      //read FLEX offset
      o1 = convertTo<double>(args.Read("o1"));

      npoly = convertTo<int>(args.ReadWithDefault("npoly","1"));

      //inference options
      inferdw = true; //args.ReadBool("inferdw");
      inferk = true; //args.ReadBool("inferk");

      ardon = args.ReadBool("ardon");

      ffton = args.ReadBool("fft");

    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    

    ncomp = compspec.Nrows();
    ntpts = tevol.Nrows();

    /*
    //ARD on baseline parameters (except DC term)
    if (npoly>1) {
    for (int i=1; i<npoly; i++) {
      ard_index.push_back(i+1);
    }
    }
    */
    if (ardon) {
    //ARD on PTR values
    for (int i=1 ; i<ncomp; i++) { //always assume we have first component (which should be water)
      ard_index.push_back(npoly+i+1);
    }
        }
}

void FLEXFwdModel::ModelUsage()
{ 
  cout << "\nUsage info for --model=FLEX:\n"
       << "Undefined\n"
    ;
}

void FLEXFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
  //cout << vec.t() << endl;
}

void FLEXFwdModel::NameParams(vector<string>& names) const
{
  names.clear();

  // name the parameters for the pools using letters
  string lettervec [] = {"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"};
  
  // baseline
  for (int i=1; i<=npoly; i++) {
    names.push_back("BL_" + lettervec[i-1]);
  }
  
  // components
    for (int i=1; i<=ncomp; i++) {
      names.push_back("PTR_" + lettervec[i-1]);
    }

for (int i=1; i<=ncomp; i++) {
      names.push_back("deltw_" + lettervec[i-1]);
    }
  
    for (int i=1; i<=ncomp; i++) {
      names.push_back("kevol_" + lettervec[i-1]);
    }

    for (int i=1; i<=ncomp; i++) {
      names.push_back("phs_" + lettervec[i-1]);
    }


}

void FLEXFwdModel::SetupARD( const MVNDist& theta, MVNDist& thetaPrior, double& Fard)
{
  Tracer_Plus tr("FLEXFwdModel::SetupARD");

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

void FLEXFwdModel::UpdateARD(
				const MVNDist& theta,
				MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("FLEXFwdModel::UpdateARD");
  
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
      //      if (b>1e12) b=1e12;
      Fard += -1.5*(log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5*log(b); //taking c as 0.5 - which it will be!
    }
  }

  return;

  }

