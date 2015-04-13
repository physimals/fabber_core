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
  return "$Id: fwdmodel_asl_multiphase.cc,v 1.3 2014/11/07 11:43:43 chappell Exp $";
}

void MultiPhaseASLFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("MultiPhaseASLFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e12;

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

     //flow vel
     if (incvel) {
       prior.means(4) = 0.3;
       if (infervel) {
	 precisions(4,4) = 10;
       }
     }

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

//init the mid phase value - by finding the point where the max intensity is
  int ind;
  float val;
  val = dmean.Maximum1(ind); // find the max
  val = (ind-1)*180/M_PI; //frequency of the minimum in ppm
  if (val>179) val -= 360;
  val *= M_PI/180;
  posterior.means(2) = val;

  if (incvel) {
    posterior.means(4) = 0.3;
  }
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
   double mag;
   double phase;
   double phaserad;
   double offset;
   double flowvel;

   mag=params(1);
   phaserad = params(2);  //in radians
   phase=params(2)*180/M_PI; //in degrees
   offset=params(3);

   if (incvel) {
     flowvel = params(4);
   }
   else {
     flowvel = 0.3;
   }   

   int nn=8*repeats;
   result.ReSize(nn);
   // loop to create result
   for (int i=1; i<=8; i++)
     {
       double evalfunc;

       double ph = (360/8*(i-1) ); //in degrees
       if (ph>179) ph -= 360; 
       double ph_rad = ph*M_PI/180; //in radians

       if (modfn == "fermi") {
	 // use the Fermi modulation function
	 evalfunc = mag*( -2/(1+exp( (abs(ph - phase) -alpha)/beta)) ) + offset; //note using the given values requires phases here to be in degrees
       }
       else if (modfn == "mat") {
	 //evaluation modulation function from interpolation of values
	 evalfunc = mag*(mod_fn(ph_rad - phaserad,flowvel)) + offset;
       }

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
      repeats = convertTo<int>(args.ReadWithDefault("repeats","1")); // number of repeats in data

      // modulation function
      modfn = args.ReadWithDefault("modfn","fermi");      

      //modmat
      string modmatstring;
      Matrix mod_temp;
      modmatstring = args.ReadWithDefault("modmat","none");
      if (modmatstring != "none") {
	mod_temp = read_ascii_matrix(modmatstring);
      }

      // shape of the fermi function
      alpha = convertTo<double>(args.ReadWithDefault("alpha","55"));
      beta = convertTo<double>(args.ReadWithDefault("beta","12"));
      
      // deal with ARD selection
      //doard=false;
      //if (inferart==true && ardoff==false) { doard=true; }

      infervel=false;
      incvel=false;
      if ( modfn == "mat" ) {
	assert(mod_temp(1,1)==99);
	int nphasepts = mod_temp.Nrows()-1;
	nvelpts = mod_temp.Ncols()-1;
	
	mod_phase = ( mod_temp.SubMatrix(2,nphasepts+1,1,1) ).AsColumn();
	mod_v = ( mod_temp.SubMatrix(1,1,2,nvelpts+1) ).AsColumn();
	mod_mat = mod_temp.SubMatrix(2,nphasepts+1,2,nvelpts+1);
	
	vmax = mod_v(nvelpts);
	vmin = mod_v(1);
	
	
	infervel = args.ReadBool("infervel");
	if (infervel) {
	  incvel=true;
	}
	else {
	  incvel = args.ReadBool("incvel");
	}
	
      }

     
      // add information about the parameters to the log
      // test correctness of specified modulation function
      if (modfn == "fermi") {
	LOG << "Inference using Fermi model" << endl;
	LOG << "alpha=" << alpha << " ,beta=" << beta << endl;
      }
      else if (modfn == "mat") {
	LOG << "Inference using numerical modulation function" << endl;
	LOG << "File is: " << modmatstring << endl;
      }
      else {
	throw invalid_argument("Unrecognised modulation function");
      }
      
	  
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
  if (incvel) {
    names.push_back("vel");
  }
}

double MultiPhaseASLFwdModel::mod_fn( const double inphase, const double v) const {
    Tracer_Plus trace("MultiPhaseASLFwdModel::mod_fn");
    // the modulation function - evaluated from interpolation of modmat
    double ans;
    double phase=inphase;

    // phase will be normally in range 0 --> 2*pi
    if (phase<0.0) phase = 0.0;
    else if (phase>2*M_PI) phase = 2*M_PI;
      // old from veasl model
      //deal with phase outside range -pi --> +pi
      //phase = asin(sin(phase)); //this assumes symmtery of function
      //if (phase>0) phase=std::fmod(phase+M_PI,2*M_PI)-M_PI;
      //else if (phase<0) phase=std::fmod(phase-M_PI,2*M_PI)+M_PI;
     // ** end old

      //bilinear interpolation
      if (v >= vmax) {
	ColumnVector usecolumn = mod_mat.Column(nvelpts);
	ans = interp(mod_phase,usecolumn,phase);
      }
      else if (v <= vmin) {
	ColumnVector usecolumn = mod_mat.Column(1);
	ans = interp(mod_phase,usecolumn,phase);
      }
      else {
	int ind=1;      
	while (v >= mod_v(ind)) ind++;

	ColumnVector usecolumn = mod_mat.Column(ind-1);
	double mod_l = interp(mod_phase,usecolumn,phase);
	ColumnVector usecolumn2 = mod_mat.Column(ind);
	double mod_u = interp(mod_phase,usecolumn2,phase);
	ans = mod_l + (v - mod_v(ind-1))/(mod_v(ind) - mod_v(ind-1)) * (mod_u-mod_l);
      }

    return ans;
  }

  double MultiPhaseASLFwdModel::interp(const ColumnVector& x, const ColumnVector& y, const double xi) const
// Look-up function for data table defined by x, y
// Returns the values yi at xi using linear interpolation
// Assumes that x is sorted in ascending order
  // ? could be replaced my MISCMATHS:interp1 ?
{
  
  double ans;
  if(xi >= x.Maximum()) 
    ans = y(x.Nrows());
  else
    if(xi <= x.Minimum()) 
      ans = y(1); 
    else{
      int ind=1;
      while(xi >= x(ind))
	ind++;      
      double xa = x(ind-1), xb = x(ind), ya = y(ind-1), yb = y(ind);
      ans = ya + (xi - xa)/(xb - xa) * (yb - ya);
    }
  return ans;
}
