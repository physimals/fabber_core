/*  fwdmodel_asl_satrecov.h - Saturation Recovery curve calibration for ASL

    Michael Chappell, IBME & FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class SatrecovFwdModel : public FwdModel {
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
  { return (LFAon?4:3);  } 

  virtual ~SatrecovFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  //virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard);
  //virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;

  // Constructor
  SatrecovFwdModel(ArgsType& args);


protected: // Constants

  // Lookup the starting indices of the parameters

  // vector indices for the parameters to expereicne ARD
  vector<int> ard_index;


  // scan parameters
  int repeats;
  int nphases;
  double t1;
  double slicedt;

  double FAnom;
  double LFA;
  double dti;
  float dg;

  bool looklocker;
  bool LFAon;
  bool fixA;

  ColumnVector tis;
  Real timax;


};
