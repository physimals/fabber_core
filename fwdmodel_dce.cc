/*  fwdmodel_dec.cc - Implements the Dynamic Contrast Enhanced MRI model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce.h"
#include <cmath>
#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

//void BuxtonFwdModel::Evaluate(volume4D<float> wholeimage<float>, volume<float> Kev, volume<float> T10, volume<float> Ve, 
//volume<float> oi3d, volume4d<float> * finalimage, volume4D<float> * defx,volume4D<float> * defy,volume4D<float> * defz)


void DCEFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("DCEFwdModel::Evaluate");

  float Kep, T10, Ve, oi3d, onset_time;
  //int tdim;
  Kep=params(1);
  T10=params(2);
  Ve=params(3);
  oi3d=params(4);
  onset_time=params(5);

if (Kep<0.001)
{ Kep= 0.001;
}

if (Kep>1)
{ Kep= 1;
}

if (oi3d<0)
{ oi3d= 0;
}

if (Ve<0.001)
{ Ve= 0.001;
}

if (Ve>1)
{ Ve= 1;
}

  // TODO - should encode the timing values somewhere (i.e. pixdim[4])
  
  // params for the matlab function model_prediction3d- provided as fixed values within the function
  //float dosage=0.1; 
  //float R1= 4.5;
  //float R2= 5.5;
  //float alpha= 3.14159265*12.0/180.0;
  //float TR= 0.0045;
  //float TE= 0.0022;

  // params for concentration_orton_3d
  //float B= 323;
  //float G= 1.07;
  //float ub= 20.2;
  //float ug= 0.172;
  //float ab= 344;
  //float ag= 1.24;
  //float t0= 0.0804;
	  
//same things for the pre-clinical data
//ideally this should be stored as a separate model and provided as an option.
// params for the matlab function model_prediction3d- provided as fixed values within the function

  float dosage=1.2; 
  float R1= 3.9;
  float R2= 4.3;
  float alpha= 3.14159265*5.0/180.0;
  float TR= 0.0011;
  float TE= 0.000585;

  // params for concentration_orton_3d
  float B= 323;
  float G= 1.07;
  float ub= 20.2;
  float ug= 0.172;
  float ab= 344;
  float ag= 1.24;
  float t0= 0.0804;


  ColumnVector conc(tdim);
  ColumnVector Q(tdim);
  ColumnVector ModelPrediction(tdim);
  ColumnVector time(tdim);
  float P;
  
  P= TR/T10;
  //time starting from zero here:
  for( int t=1; t<=tdim; t++) {
    time(t) = 0.29*t-0.29*onset_time-0.29;
// cout<<time<<" "<<endl;

	if (time(t)<0.001)
	{ time(t)= 0.001;
	}
    
    // Orton AIF for mouse
     conc(t)= dosage*Kep*Ve*((5.8*(exp(-Kep*time(t))-exp(-1.78*time(t)))/(1.78-Kep)) +
     			    (3.1*(exp(-Kep*time(t))-exp(-0.045*time(t)))/(0.045-Kep)));
    
    // ORTON AIF
    // conc(t)= B*Kep*Ve*(time(t)*exp(-ub*time(t))+(exp(-Kep*time(t))-exp(-ub*time(t)))/(Kep-ub))/(Kep-ub)+
      //  G*Kep*Ve*((exp(-ug*time(t))-exp(-Kep*time(t)))/(Kep-ug)-(exp(-ub*time(t))-exp(-Kep*time(t)))/(Kep-ub));
    
    Q(t)= R1*TR*conc(t);
    
    ModelPrediction(t) = oi3d*exp(-R2*TE*conc(t))*(1-exp(-P-Q(t))-cos(alpha)*
						   (exp(-P)-exp(-2*P-Q(t))))/(1-exp(-P)-cos(alpha)*
									      (exp(-P-Q(t))-exp(-2*P-Q(t))));
//cout<<ModelPrediction(t)<<" "<<endl;
    
    // This final output should be a vector of the same size as 'time' 
  }
	
  result=ModelPrediction;

}

string DCEFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_dce.cc,v 1.5 2012/11/29 16:44:47 chappell Exp $";
}

void DCEFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("DCEFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());


    // set prior, mean and precision for each parameter

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // Kep
     prior.means(1) = 0.5;
     precisions(1,1) = 0.5;

    // T10
     prior.means(2) = 1;
     precisions(2,2) = 1000000;

    // Ve
     prior.means(3) = 0.5;
     precisions(3,3) = 1;

    // intensity of first volume
     prior.means(4) = 0.15;
     precisions(4,4) = 10000;

    // onset time
     prior.means(5) = 11;
     precisions(5,5) = 5;
    
    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;

    // Modify posterior for perfusion to a more realisitic starting point
    // set initial guess for uninformative priors
    posterior.means(4) = 0.1;
    precisions(4,4) = 100;
    posterior.SetPrecisions(precisions);
    
}    

DCEFwdModel::DCEFwdModel(ArgsType& args)
{
    string scanParams = args.ReadWithDefault("scan-params","cmdline");


    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
      tdim = convertTo<int>(args.Read("tdim")); // number of repeats in data

      // add information about the parameters to the log
      // LOG << "    Data parameters: #repeats = " << repeats << ", t1 = " << t1 << ", t1b = " << t1b;
	  
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    
    
    //turn off all ARD for timebeing
    doard=false;
}

void DCEFwdModel::ModelUsage()
{ 
  cout << "\nUsage info for --model=dce:\n"  // TODO - need to set this up at a higher level
       << "Required parameters:\n"
       << "--tdim=<no. time points>\n"
    ;
}

void DCEFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void DCEFwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  names.push_back("Kep");
  names.push_back("T10");
  names.push_back("Ve");
  names.push_back("first_vol_intensity");
  names.push_back("onset_time");
}

// LEAVE THIS AS IS
void DCEFwdModel::SetupARD( const MVNDist& theta, MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("DCEFwdModel::SetupARD");

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

// LEAVE AS IS
void DCEFwdModel::UpdateARD(
				const MVNDist& theta,
				MVNDist& thetaPrior, double& Fard) const
{
  Tracer_Plus tr("DCEFwdModel::UpdateARD");
  
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



