/*  noisemodel_ar.cc - Class implementation for the AR(1) noise model

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

    Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */
 
#include "noisemodel.h"
#include "dist_gamma.h"
#include <vector>

using namespace std;

class Ar1cParams;

// Helper class -- caches some of the AR matrices
class Ar1cMatrixCache {
public:
  const SymmetricBandMatrix& GetMatrix(unsigned n, unsigned a12pow, 
                       unsigned a3pow) const;
  const SymmetricBandMatrix& GetMarginal(unsigned n) const;

  void Update(const Ar1cParams& dist, int nTimes);

  Ar1cMatrixCache(int numPhis) : nPhis(numPhis) { return; }
  Ar1cMatrixCache(const Ar1cMatrixCache& from)
    : alphaMarginals(from.alphaMarginals), alphaMatrices(from.alphaMatrices),
      nPhis(from.nPhis)
    { return; }

private:
  vector<SymmetricBandMatrix> alphaMarginals; 
       // recalculated whenever alpha changes
  unsigned FlattenIndex(unsigned n, unsigned a12pow, unsigned a34pow) const
    { assert(n==1 || n==2 && a12pow<=2 && a34pow<=2);
      return n-1 + 2*( a12pow + 3*(a34pow) ); } 

  vector<SymmetricBandMatrix> alphaMatrices; 
       // should only be calculated once
       // Note that if more than one model is being inferred upon at a time,
       // this will be unnecessarily duplicated in every one of them --
       // might speed things up considerably by sharing.
       
  int nPhis;
};

 
// Parameter-storage class -- it's really just an enhanced structure
class Ar1cParams : public NoiseParams {
public:
    virtual Ar1cParams* Clone() const
        { return new Ar1cParams(*this); } 

    virtual const Ar1cParams& operator=(const NoiseParams& in)
      { const Ar1cParams& from = dynamic_cast<const Ar1cParams&>(in);
	alpha = from.alpha; phis = from.phis; alphaMat = from.alphaMat; return *this; }

    virtual const MVNDist OutputAsMVN() const;
    virtual void InputFromMVN(const MVNDist& mvn);
       
    // Human-readable debug output (dump internal state to LOG)
    virtual void Dump(const string indent = "") const;

    // Constructor/destructor
    Ar1cParams(int nAlpha, int nPhi) : 
        alpha(nAlpha), phis(nPhi), alphaMat(nPhi) { return; }
    Ar1cParams(const Ar1cParams& from) : 
        alpha(from.alpha), phis(from.phis), alphaMat(from.alphaMat) { return; }
    virtual ~Ar1cParams() { return; }

private:
    friend class Ar1cNoiseModel; // Needs to use this class like it's a structure
    friend class Ar1cMatrixCache;
    MVNDist alpha;
    vector<GammaDist> phis;
    
    Ar1cMatrixCache alphaMat;
};


class Ar1cNoiseModel : public NoiseModel {
 public:
  static NoiseModel* NewInstance();

//  virtual Ar1cNoiseModel* Clone() const;
  // makes a new identical copy of this object
  
    virtual Ar1cParams* NewParams() const
        { return new Ar1cParams( NumAlphas(), nPhis ); }

    virtual void HardcodedInitialDists(NoiseParams& prior, NoiseParams& posterior) const;


//  virtual void LoadPrior( const string& filename );
  // loads priors from file, and also initializes posteriors
  
    virtual void Precalculate( NoiseParams& noise, const NoiseParams& noisePrior,
        const ColumnVector& sampleData ) const;
    // Used to pre-evaluate the alpha matrices in the cache

  // virtual void AdjustPrior(...) might be needed for multi-voxel methods...
  // probably best for that to go in a derived class. 
  
//  virtual void Dump(const string indent = "") const;
//  virtual void DumpPrior(const string indent = "") const;  
//  virtual void DumpPosterior(const string indent = "") const; 
  // human-readable debug output

//  virtual const MVNDist GetResultsAsMVN() const;

  // Constructor/destructor
    virtual void Initialize(FabberRunData& args);
    //Ar1cNoiseModel(const string& ar1CrossTerms, int numPhis );
    // ar1CrossTerms must be either "none", "dual", or "same". 
    
    virtual ~Ar1cNoiseModel() { return; }
  
  // VB Updates
    
  virtual void UpdateNoise(
    NoiseParams& noise, 
    const NoiseParams& noisePrior, 
  	const MVNDist& theta,
  	const LinearFwdModel& linear,
  	const ColumnVector& data) const
  { UpdateAlpha(noise, noisePrior, theta, linear, data);
    UpdatePhi(noise, noisePrior, theta, linear, data); } 

  virtual void UpdateAlpha(
    NoiseParams& noise, 
    const NoiseParams& noisePrior, 
    const MVNDist& theta,
    const LinearFwdModel& model,
    const ColumnVector& data) const;
    
  virtual void UpdatePhi(
    NoiseParams& noise, 
    const NoiseParams& noisePrior,   
    const MVNDist& theta,
    const LinearFwdModel& model,
    const ColumnVector& data) const;

  virtual void UpdateTheta(
    const NoiseParams& noise, 
//    const NoiseParams& noisePrior,   
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

   int NumParams();

//  void SaveParams(const MVNDist& theta) {};		 
//  void RevertParams(MVNDist& theta) {};

 protected: 
//  Ar1cParameters* prior;
//  Ar1cParameters* posterior;  
    // Whenever this changes, call alphaMat.Update!

//  Ar1cMatrixCache alphaMat;
  string ar1Type;
  int NumAlphas() const; // converts the above string into a number
  int nPhis;
};

