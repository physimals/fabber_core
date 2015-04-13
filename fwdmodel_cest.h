/*  fwdmodel_cest_devel.h - Development resting state ASL model

    Michael Chappell, IBME PUMMA & FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class CESTFwdModel : public FwdModel {
public: 
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;
  void Initialise(MVNDist& posterior) const;

   static void ModelUsage();
  virtual string ModelVersion() const;
                  
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const 
  { return (3*npool) + 1 + ( inferdrift?1:0 ) + ( t12soft? (2*npool):0 );
  } 

  virtual ~CESTFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  using FwdModel::SetupARD;
  virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard);
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;


  // Constructor
  CESTFwdModel(ArgsType& args);


protected: 
  //specific functions
  void Mz_spectrum(ColumnVector& Mz, const ColumnVector& wvec, const ColumnVector& w1, const ColumnVector& t, const ColumnVector& M0, const Matrix& wi, const Matrix& kij, const Matrix& T12) const;
  ReturnMatrix Mz_spectrum_lorentz(const ColumnVector& wvec, const ColumnVector& w1, const ColumnVector& t, const ColumnVector& M0, const Matrix& wi, const Matrix& kij, const Matrix& T12) const;

  //maths functions
  void Ainverse(const Matrix A, RowVector& Ai) const;
  ReturnMatrix expm(Matrix inmatrix) const;
  ReturnMatrix expm_eig(Matrix inmatrix) const;
  ReturnMatrix expm_pade(Matrix inmatrix) const;
  ReturnMatrix PadeApproximant(Matrix inmatrix, int m) const;
  ReturnMatrix PadeCoeffs(int m) const;

// Constants

  // Lookup the starting indices of the parameters
  

  // vector indices for the parameters to expereicne ARD
  vector<int> ard_index;

  //flags
  bool t12soft;
  bool inferdrift;
  //bool pvcorr;
  bool lorentz;
  bool steadystate;

  // scan parameters
  vector<float> t;

  //model parameters
  int npool;
  float wlam;
  float B1set;
  ColumnVector poolppm;
  ColumnVector poolk;
  Matrix T12master;
  //Matrix T12WMmaster;
  //Matrix T12CSFmaster;

  // Data specification
  ColumnVector wvec;
  ColumnVector w1vec;
  ColumnVector tsatvec;

  //pulse specificaiton
  ColumnVector pmagvec;
  ColumnVector ptvec;
  int nseg;

  // processing flags
  mutable bool fastgrad; //use a fast approximation to the expm because we are caculating the gradient
  
  // ard flags
  bool doard;
 
};
