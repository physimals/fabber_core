/*  noisemodel_white.cc - Class implementation for the multiple white noise model

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "noisemodel_white.h"
#include "noisemodel.h"
#include <stdexcept>
#include "miscmaths/miscmaths.h"
using namespace MISCMATHS;
using namespace Utilities;
#include "easylog.h"

void WhiteNoiseModel::HardcodedInitialDists(NoiseParams& priorIn,
    NoiseParams& posteriorIn) const
{
    Tracer_Plus tr("WhiteNoiseModel::HardcodedInitialDists");
    
    WhiteParams& prior = dynamic_cast<WhiteParams&>(priorIn);
    WhiteParams& posterior = dynamic_cast<WhiteParams&>(posteriorIn);
    
    int nPhis = Qis.size();
    assert(nPhis > 0);
//    prior.resize(nPhis);
//    posterior.resize(nPhis);
    assert(nPhis == prior.nPhis);
    assert(nPhis == posterior.nPhis);

    for (int i = 1; i <= nPhis; i++)
      {
	if (phiprior == -1) {
	  // use non-informative prior on phi
	   posterior.phis[i-1].b = prior.phis[i-1].b = 1e6;
	   posterior.phis[i-1].c = prior.phis[i-1].c = 1e-6;

       // Bit of black magic... tiny initial noise precision seems to help
	   posterior.phis[i-1].b = 1e-8;
	   posterior.phis[i-1].c = 50;
	}
	else {
	  // use input nosie std dev to determine the prior (and inital values) for phi
	  // this is for the case where N is small, so it make sense to use a large(ish) c_prior
	  // since C ~ N (data points) and we struggle to estimate noise when N is small.
	  posterior.phis[i-1].c = prior.phis[i-1].c = 0.5; //assumes that a given std deviation is equivelent to N=1 measurements
	  posterior.phis[i-1].b = prior.phis[i-1].b = 1/(phiprior*phiprior*prior.phis[i-1].c); //NB phiprior is a std dev for the noise that is read in
	}
    }
}

const MVNDist WhiteParams::OutputAsMVN() const
{
  Tracer_Plus tr("WhiteParams::OutputAsMVN");
  
  assert((unsigned)nPhis == phis.size());
  MVNDist mvn( phis.size() );
  SymmetricMatrix vars( phis.size() );
  vars = 0;
  for (unsigned i = 1; i <= phis.size(); i++)
    {
      mvn.means(i) = phis[i-1].CalcMean();
      vars(i,i) = phis[i-1].CalcVariance();
    }
  mvn.SetCovariance(vars);
  return mvn;
}

void WhiteParams::InputFromMVN( const MVNDist& mvn )
{
    for (unsigned i = 1; i <= phis.size(); i++)
    {
        phis[i-1].SetMeanVariance( mvn.means(i), mvn.GetCovariance()(i,i) );
        for (int j = i+1; j <= mvn.means.Nrows(); j++)
            if (mvn.GetCovariance()(i, j) != 0.0)
                throw Invalid_option("Phis should have zero covariance!");
    }    
} 


void WhiteParams::Dump(const string indent) const
{
  Tracer_Plus tr("WhiteParams::Dump");
  assert( (unsigned)nPhis == phis.size() );
  for (unsigned i = 0; i < phis.size(); i++)
    {
      LOG << indent << "Phi_" << i+1 << ": ";
      phis[i].Dump();
    }
}

//WhiteNoiseModel::WhiteNoiseModel(const string& pattern) 
//  : phiPattern(pattern)
WhiteNoiseModel::WhiteNoiseModel(ArgsType& args)
  : phiPattern(args.ReadWithDefault("noise-pattern","1"))
{ 
  Tracer_Plus tr("WhiteNoiseModel::WhiteNoiseModel");
  assert(phiPattern.length() > 0);
  MakeQis(phiPattern.length()); // a quick way to validate the input

  // Allow phi to be locked externally
  lockedNoiseStdev = convertTo<double>(args.ReadWithDefault("locked-noise-stdev","-1"));
  assert(lockedNoiseStdev == -1 || lockedNoiseStdev > 0); // TODO (should throw)

  // Set phi prior externally
  phiprior = convertTo<double>(args.ReadWithDefault("prior-noise-stddev","-1"));
  if ( phiprior < 0 && phiprior != -1 ) throw Invalid_option("Incorrect choice of prior-noise-stddev, is should be >0");
}


void WhiteNoiseModel::MakeQis(int dataLen) const
{
  Tracer_Plus tr("WhiteNoiseModel::MakeQis");
  if (!Qis.empty() && Qis[0].Nrows() == dataLen) 
    return;  // Qis are already up-to-date

  // Read the pattern string into a vector pat
  const int patternLen = phiPattern.length();
  assert(patternLen > 0);
  if (patternLen > dataLen)
    throw Invalid_option(
      "Pattern length exceeds data length... this is probably a mistake");

  vector<int> pat;
  int nPhis = 0;

  for (int i = 1; i <= patternLen; i++)
    {
      char c = phiPattern[i-1];
      int n;
      if (c >= '1' && c <= '9') // Index from 1, so 0 is not allowed
	n = (c - '0');
      else if (c >= 'A' && c <= 'Z')
	n = (c - 'A' + 10);
      else if (c >= 'a' && c <= 'z')
	n = (c - 'a' + 10);
      else
	throw Invalid_option(
	          string("Invalid character in pattern: '") + c + "'");
      pat.push_back(n);
      if (nPhis<n) nPhis = n;
    }


  //  cout << "Pat is " << pat << endl;

  // Extend pat to the the full data length by repeating
  while ((int)pat.size() < dataLen)
    pat.push_back(pat.at(pat.size()-patternLen));

  LOG << "Pattern of phis used is " << pat << endl;

  // Regenerate the Qis
  Qis.clear();
  DiagonalMatrix zeroes(dataLen);
  zeroes = 0.0;
  Qis.resize(nPhis, zeroes); // initialize all to zero

  for (int d = 1; d <= dataLen; d++)
    {
    //    Qis[pat[d-1]-1](d,d) = 1;
    Qis.at(pat.at(d-1)-1)(d,d) = 1;
    }

  // Sanity checking
  for (int i = 1; i <= nPhis; i++)
    if (Qis[i-1].Trace() < 1.0) // this phi is never used
      throw Invalid_option(
	 "At least one Phi was unused! This is probably a bad thing.");
}

void WhiteNoiseModel::UpdateNoise(
    NoiseParams& noise,
    const NoiseParams& noisePrior,
  	const MVNDist& theta,
  	const LinearFwdModel& linear,
  	const ColumnVector& data) const
{
  Tracer_Plus tr("WhiteNoiseModel::UpdateNoise");
  
  WhiteParams& posterior = dynamic_cast<WhiteParams&>(noise);
  const WhiteParams& prior = dynamic_cast<const WhiteParams&>(noisePrior);
  
  const Matrix& J = linear.Jacobian();
  ColumnVector k = data - linear.Offset() + J*(linear.Centre() - theta.means);
  
  // check the Qis are valid
  MakeQis(data.Nrows());
  const int nPhis = Qis.size();
  assert(nPhis == posterior.nPhis);
  assert(nPhis == prior.nPhis);

  // Update each phi distribution in turn
  for (int i = 1; i <= nPhis; i++)
    {
      const DiagonalMatrix& Qi = Qis[i-1];
      double tmp = 
	(k.t() * Qi * k).AsScalar()
	+ (theta.GetCovariance() * J.t() * Qi * J).Trace();
      
      posterior.phis[i-1].b =
	1/( tmp*0.5 + 1/prior.phis[i-1].b);
      
      double nTimes = Qi.Trace(); // number of data points for this dist
      assert(nTimes == int(nTimes)); // should be an integer

      posterior.phis[i-1].c = 
	(nTimes-1)*0.5 + prior.phis[i-1].c;

      if (lockedNoiseStdev > 0)
	{
	  // Ignore this update and force phi to a specified value.
	  // b*c = noise precision = lockedNoiseStdev^-2
	  posterior.phis[i-1].b = 
	    1/posterior.phis[i-1].c/lockedNoiseStdev/lockedNoiseStdev;
	}

    }
}

void WhiteNoiseModel::UpdateTheta(
    const NoiseParams& noiseIn,
  	MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& linear,
        const ColumnVector& data,
        MVNDist* thetaWithoutPrior,
	float LMalpha) const
{
  Tracer_Plus tr("WhiteNoiseModel::UpdateTheta");

  //cout << "start:" << theta.means.t() << endl;

  //  if (thetaWithoutPrior != NULL)
  //    throw Invalid_option("This noise model doesn't yet work with "
  //      +string("--spatial-prior-output-correction... implement it!\n"));

  const WhiteParams& noise = dynamic_cast<const WhiteParams&>(noiseIn);

  const ColumnVector &ml = linear.Centre();
  const ColumnVector &gml = linear.Offset();
  const Matrix &J = linear.Jacobian();

  // Make sure Qis are up-to-date
  MakeQis(data.Nrows());
  assert(Qis.size() == (unsigned)noise.nPhis);

  // Marginalize over phi distributions
  DiagonalMatrix X(data.Nrows()); 
  X = 0;
  for (unsigned i = 1; i <= Qis.size(); i++)
    X += Qis[i-1] * noise.phis[i-1].CalcMean();
  // Qis are diagonal matrices with 1 only where that phi applies.
  // Adding up all the Qis will give you the identity matrix.

  // Calculate Lambda & Lambda*m (without priors)
  SymmetricMatrix Ltmp;
  Ltmp << J.t() * X * J;  
  // use << instead of = because this is considered a lossy assignment
  // (since NEWMAT isn't smart enough to know J'*X*J is always symmetric)
  ColumnVector mTmp = J.t() * X * (data - gml + J*ml);

  // Update Lambda and m (including priors)
  theta.SetPrecisions(thetaPrior.GetPrecisions() + Ltmp);

// Error checking
  LogAndSign chk = theta.GetPrecisions().LogDeterminant();
  if (chk.Sign() <= 0) {
    LOG << "Note: In UpdateTheta, theta precisions aren't positive-definite: "
	<< chk.Sign() << ", " << chk.LogValue() << endl;
  }      

  if (LMalpha>0.0) {
    //we are in LM mode so use the appropriate update
    //cout << LMalpha << endl;
    Matrix Delta;
    SymmetricMatrix prec;
    DiagonalMatrix precdiag;
    prec = theta.GetPrecisions();
    precdiag << prec;

    // a different (but equivalent?) form for the LM update
    Delta = J.t() * X * (data - gml) + thetaPrior.GetPrecisions()*thetaPrior.means - thetaPrior.GetPrecisions()*ml;
    theta.means = ml + (prec + LMalpha*precdiag).i()*Delta;

    // LM update
    //theta.means = (prec + LMalpha*precdiag).i()
    // * ( mTmp + thetaPrior.GetPrecisions() * thetaPrior.means );

  }
  else { //normal update (NB the LM update resuces to this when alpha=0 strictly)
  theta.means = theta.GetCovariance()
    * ( mTmp + thetaPrior.GetPrecisions() * thetaPrior.means );
  }

  if (thetaWithoutPrior != NULL)
    {
      Tracer_Plus tr("WhiteNoiseModel::UpdateTheta - WithoutPrior calcs");
      thetaWithoutPrior->SetSize(theta.GetSize());
      
      thetaWithoutPrior->SetPrecisions(Ltmp);

      //      try
      //	{
      //          thetaWithoutPrior->GetCovariance(); // Cached, so not wasteful.
      //	}
      //      catch (Exception)
      //	{
      //	  Warning::IssueAlways("Ltmp was singular -- adding 1e-20 to diagonal");
      //	  thetaWithoutPrior->SetPrecisions(Ltmp + IdentityMatrix(Ltmp.Nrows()) * 1e-20);
      //	}

      thetaWithoutPrior->means = thetaWithoutPrior->GetCovariance() * mTmp;
    }

  /*
  // Error checking
  LogAndSign chk = theta.GetPrecisions().LogDeterminant();
  if (chk.Sign() <= 0)
    LOG << "Note: In UpdateTheta, theta precisions aren't positive-definite: "
	<< chk.Sign() << ", " << chk.LogValue() << endl;
  */

  //cout << "end:" << theta.means.t() << endl;
}


double WhiteNoiseModel::CalcFreeEnergy(
    const NoiseParams& noiseIn,
    const NoiseParams& noisePriorIn,
	const MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& linear,
  	const ColumnVector& data) const
{
    Tracer_Plus tr("WhiteNoiseModel::CalcFreeEnergy");
    const int nPhis = Qis.size();
    const WhiteParams& noise = dynamic_cast<const WhiteParams&>(noiseIn);
    const WhiteParams& noisePrior = dynamic_cast<const WhiteParams&>(noisePriorIn);

  // Calculate some matrices we will need
  const Matrix &J = linear.Jacobian();
  ColumnVector k = data - linear.Offset()
    + J * (linear.Centre() - theta.means);
  const SymmetricMatrix& Linv = theta.GetCovariance();

  // some values we will need
  int nTimes = data.Nrows(); //*NB assume that each row is an individual time point
  int nTheta = theta.means.Nrows();

  // The following is based on noisemodel_ar::CalcFreeEnergy, modified to remove ar parts - MAC 11-7-2007
  // Some modifications have been made for consistency with (MAC)varbayes2.m - these are noted

  // calcualte individual aprts of the free energy
  double expectedLogThetaDist = //bits arising from the factorised posterior for theta
    +0.5 * theta.GetPrecisions().LogDeterminant().LogValue()
    -0.5 * nTheta * (log(2*M_PI) + 1);

  double expectedLogPhiDist = 0; //bits arising fromt he factorised posterior for phi
  vector<double> expectedLogPosteriorParts(10); //bits arising from the likelihood
   for (int i = 0; i<10; i++) 
     expectedLogPosteriorParts[i] = 0;
  
   for (int i = 0; i < nPhis; i++) 
     {
      double si = noise.phis[i].b;
      double ci = noise.phis[i].c;
      double siPrior = noisePrior.phis[i].b;
      double ciPrior = noisePrior.phis[i].c;

      expectedLogPhiDist +=
	-gammaln(ci) - ci*log(si) - ci
	+(ci-1)*(digamma(ci)+log(si));
      
      expectedLogPosteriorParts[0] += 
	(digamma(ci)+log(si)) * ( (Qis[i].Trace())*0.5 + ciPrior - 1); // nTimes using phi_{i+1} = Qis[i].Trace()
      
      expectedLogPosteriorParts[9] += 
	-gammaln(ciPrior) -ciPrior*log(siPrior) - si*ci/siPrior;
     } 

   expectedLogPosteriorParts[1] = 0; //*NB not required
  
  expectedLogPosteriorParts[2] =
    -0.5 * (k.t() *  k).AsScalar() 
    -0.5 * (J.t() * J * Linv).Trace(); //*NB remove Qsum
  
  expectedLogPosteriorParts[3] =
    +0.5 * thetaPrior.GetPrecisions().LogDeterminant().LogValue()
    -0.5 * nTimes * log(2*M_PI)
    -0.5 * nTheta * log(2*M_PI);
  
  expectedLogPosteriorParts[4] = 
    -0.5 * (
            (theta.means - thetaPrior.means).t() 
            * thetaPrior.GetPrecisions() 
            * (theta.means - thetaPrior.means)
	    ).AsScalar();
  
  expectedLogPosteriorParts[5] =
    -0.5 * (Linv * thetaPrior.GetPrecisions()).Trace();

  expectedLogPosteriorParts[6] = 0; //*NB not required
  
  expectedLogPosteriorParts[7] = 0; //*NB not required  
  
  expectedLogPosteriorParts[8] = 0; //*NB not required

  
  // Assemble the parts into F
  double F =  
    - expectedLogThetaDist 
    - expectedLogPhiDist;
  
  for (int i=0; i<10; i++) 
    F += expectedLogPosteriorParts[i];

  // Error checking
  if (! (F - F == 0) )
  {
    LOG_ERR("expectedLogThetaDist == " << expectedLogThetaDist << endl);
    LOG_ERR("expectedLogPhiDist == " << expectedLogPhiDist << endl);
    //LOG_ERR("expectedLogPosteriorParts == " << expectedLogPosteriorParts << endl);
    throw overflow_error("Non-finite free energy!");
  }
  
  return F;
}

/*
void WhiteNoiseModel::SaveParams(const MVNDist& theta) {
  Tracer_Plus tr("WhiteNoiseModel::SaveParams");
  // save the current values of parameters 
  int nPhis = phis.size();
  assert(nPhis > 0);
  phissave.resize(phis.size());
  for (int i = 1; i <= nPhis; i++)
      {
	phissave[i-1].b = phis[i-1].b ;
	phissave[i-1].c = phis[i-1].c ;
      }
  thetasave = theta;
}

void WhiteNoiseModel::RevertParams(MVNDist& theta) {
  Tracer_Plus tr("WhiteNoiseModel::RevertParams");
  int nPhis = phis.size();
  for (int i = 1; i <= nPhis; i++)
      {
	phis[i-1].b = phissave[i-1].b ;
	phis[i-1].c = phissave[i-1].c ;
      }
  theta = thetasave;
}
*/


