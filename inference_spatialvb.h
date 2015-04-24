/*  inference_spatialvb.h - implementation of VB with spatial priors

    Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

    Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_vb.h"
#ifndef __FABBER_LIBRARYONLY
#endif //__FABBER_LIBRARYONLY

class CovarianceCache {
 public:
#ifndef __FABBER_LIBRARYONLY
  void CalcDistances(const NEWIMAGE::volume<float>& mask, const string& distanceMeasure);
#endif //__FABBER_LIBRARYONLY
  void CalcDistances(const NEWMAT::Matrix& voxelCoords, const string& distanceMeasure);
  const SymmetricMatrix& GetDistances() const { return distances; }

  const ReturnMatrix GetC(double delta) const; // quick to calculate
  const SymmetricMatrix& GetCinv(double delta) const;

  //  const Matrix& GetCiCodist(double delta) const;
  const SymmetricMatrix& GetCiCodistCi(double delta, double* CiCodistTrace = NULL) const;

  bool GetCachedInRange(double *guess, double lower, double upper, bool allowEndpoints = false) const;
  // If there's a cached value in (lower, upper), set *guess = value and 
  // return true; otherwise return false and don't change *guess.

 private:
  SymmetricMatrix distances;
  typedef map<double, SymmetricMatrix> Cinv_cache_type;
  mutable Cinv_cache_type Cinv_cache; 
  
  typedef map<double, pair<SymmetricMatrix,double> > CiCodistCi_cache_type;
  //  mutable CiCodist_cache_type CiCodist_cache; // only really use the Trace
  mutable CiCodistCi_cache_type CiCodistCi_cache;
};




class SpatialVariationalBayes : public VariationalBayesInferenceTechnique {
public:
  static InferenceTechnique* NewInstance();
    SpatialVariationalBayes() : 
        VariationalBayesInferenceTechnique(), 
        spatialDims(-1) { return; }
    virtual void Initialize(FwdModel* fwd_model, ArgsType& args);
    virtual void DoCalculations(const DataSet& data);
//    virtual ~SpatialVariationalBayes();

 protected:

    int spatialDims; // 0 = no spatial norm; 2 = slice only; 3 = volume
    bool continuingFromFile;

    //    bool useDataDrivenSmoothness;
    //    bool useShrinkageMethod;
    //    bool useDirichletBC;
    //    bool useMRF;
    //    bool useMRF2; // without the dirichlet bcs

    double maxPrecisionIncreasePerIteration; // Should be >1, or -1 = unlimited

    vector<vector<int> > neighbours; // Sparse matrix would be easier
    vector<vector<int> > neighbours2; // Sparse matrix would be easier
#ifndef __FABBER_LIBRARYONLY
    void CalcNeighbours(const NEWIMAGE::volume<float>& mask);
#endif //__FABBER_LIBRARYONLY
    void CalcNeighbours(const Matrix& voxelCoords);

    //vector<string> imagepriorstr; now inherited from spatialvb
    
    // For the new (Sahani-based) smoothing method:    
    CovarianceCache covar;
    string distanceMeasure;

    double fixedDelta;
    double fixedRho;
    bool updateSpatialPriorOnFirstIteration;  
    bool useEvidenceOptimization;
    double alwaysInitialDeltaGuess;
    
    bool useFullEvidenceOptimization;
    bool useSimultaneousEvidenceOptimization;
    int firstParameterForFullEO;
    bool useCovarianceMarginalsRatherThanPrecisions;
    bool keepInterparameterCovariances;

    int newDeltaEvaluations;

    string spatialPriorsTypes; // one character per parameter
    //    bool spatialPriorOutputCorrection;

    bool bruteForceDeltaSearch;

    double OptimizeSmoothingScale(
      const DiagonalMatrix& covRatio,
      //const SymmetricMatrix& covRatioSupplemented,
      const ColumnVector& meanDiffRatio, 
      double guess, double* optimizedRho = NULL, 
      bool allowRhoToVary = true,
      bool allowDeltaToVary = true) const;

    double OptimizeEvidence(
      // const vector<MVNDist>& fwdPriorVox, // used for parameters other than k
      const vector<MVNDist*>& fwdPosteriorWithoutPrior, // used for parameter k
      // const vector<SymmetricMatrix>& Si,
      int k, const MVNDist* initialFwdPrior, double guess,
      bool allowRhoToVary = false,
      double* rhoOut = NULL) const;
};


#ifdef __FABBER_LIBRARYONLY_TESTWITHNEWIMAGE
#include "newimage/newimage.h"
using namespace NEWIMAGE;
// Helper function, useful elsewhere:
void ConvertMaskToVoxelCoordinates(const volume<float>& mask, Matrix& voxelCoords);
#endif //__FABBER_LIBRARYONLY_TESTWITHNEWIMAGE
