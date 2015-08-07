/*  noisemodel_ar.h - Class declaration for the multiple white noise model

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "noisemodel.h"
#include "dist_gamma.h"

class WhiteParams : public NoiseParams {
public:
    virtual WhiteParams* Clone() const
        { return new WhiteParams(*this); }
    virtual const WhiteParams& operator=(const NoiseParams& in)
      { const WhiteParams& from = dynamic_cast<const WhiteParams&>(in);
	assert(nPhis == from.nPhis); phis = from.phis; return *this; }
    
    virtual const MVNDist OutputAsMVN() const;
    virtual void InputFromMVN(const MVNDist& mvn);
    
    virtual void Dump(const string indent = "") const;
    
    WhiteParams(int N) : nPhis(N), phis(N) { return; }
    WhiteParams(const WhiteParams& from) 
        : nPhis(from.nPhis), phis(from.phis) { return; }
    
private:
    friend class WhiteNoiseModel;
    const int nPhis;
    vector<GammaDist> phis;
};    

class WhiteNoiseModel : public NoiseModel {
 public:

    virtual WhiteParams* NewParams() const
        { return new WhiteParams( Qis.size() ); }

    virtual void HardcodedInitialDists(NoiseParams& prior, 
        NoiseParams& posterior) const; 


  // Constructor/destructor
    //  WhiteNoiseModel(const string& pattern);
    WhiteNoiseModel(ArgsType& args);
  // pattern says which phi distribution to use for each data points; this
  // string is repeated as necessary to make up the data length. e.g. for 
  // dual-echo data, pattern = "12".  after 9, use letters (A=10, ...)
  // Simplest case: pattern = "1".

  virtual ~WhiteNoiseModel() { return; }

  // Do all the calculations
  virtual void UpdateNoise(
    NoiseParams& noise,
    const NoiseParams& noisePrior,  
  	const MVNDist& theta,
  	const LinearFwdModel& model,
  	const ColumnVector& data) const;

  virtual void UpdateTheta(
    const NoiseParams& noise,
  	MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& model,
        const ColumnVector& data,
    MVNDist* thetaWithoutPrior = NULL,
        float LMalpha = 0
    ) const;

  virtual double CalcFreeEnergy(
    const NoiseParams& noise,
    const NoiseParams& noisePrior,
	const MVNDist& theta,
  	const MVNDist& thetaPrior,
  	const LinearFwdModel& model,
  	const ColumnVector& data) const;

 protected:
  const string phiPattern;

  double lockedNoiseStdev; // Allow phi to be locked externally
  double phiprior; //allow external setting of the prior nosie std deviation (and thence phi)

  // Diagonal matrices, indicating which data points use each phi
  mutable vector<DiagonalMatrix> Qis; // mutable because it's used as a cache
  void MakeQis(int dataLen) const;
};
