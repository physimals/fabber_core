/*  fwdmodel_cest_devel.h - Development resting state ASL model

    Michael Chappell, IBME PUMMA & FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class CESTDevelFwdModel : public FwdModel {
public: 
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;
  static void ModelUsage();
  virtual string ModelVersion() const;
                  
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const 
  { return 9 + (pvcorr?9:0) + (t12soft?6:0) + (wassron?1:0);
  } 

  virtual ~CESTDevelFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard);
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;


  // Constructor
  CESTDevelFwdModel(ArgsType& args);


protected: 
  //specific functions
  ReturnMatrix Mz_spectrum(ColumnVector wvec, float w1, float t, ColumnVector M0, ColumnVector wi, Matrix kij, Matrix T12) const;
  //ReturnMatrix expm(Matrix inmatrix) const;
  //ReturnMatrix PadeApproximant(Matrix inmatrix, int m) const;
  //ReturnMatrix PadeCoeffs(int m) const;

// Constants

  // Lookup the starting indices of the parameters
  

  // vector indices for the parameters to expereicne ARD
  vector<int> ard_index;

  //flags
  bool wassron;
  bool wassronly;
  bool t12soft;
  bool pvcorr;
  bool basic;
  bool mton;

  // scan parameters
  vector<float> t;

  //model parameters
  int npool;
  vector<float> wlam;
  float W_wlam;
  float ppm_apt_set;
  //float ppm_mt;
  vector<float> B1set;
  Matrix T12master;
  Matrix T12WMmaster;
  Matrix T12CSFmaster;
  vector<ColumnVector> wvec;

  //WASSR
  float W_w1;
  float W_t;
  ColumnVector W_wvec;
  
  // ard flags
  bool doard;
 
};
