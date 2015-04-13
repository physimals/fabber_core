/*  fwdmodel_flex.h - FLEX model

    Michael Chappell, IBME PUMMA & FMRIB Image Analysis Group

    Copyright (C) 2011 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class FLEXFwdModel : public FwdModel {
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
  { return npoly + 4*ncomp ;
  } 

  virtual ~FLEXFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard);
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;


  // Constructor
  FLEXFwdModel(ArgsType& args);


protected: 
  //specific functions

// Constants

  // Lookup the starting indices of the parameters
  

  // vector indices for the parameters to expereicne ARD
  vector<int> ard_index;

  //flags

  // scan parameters
  ColumnVector tevol;

  //model parameters
  int ncomp; //number of components
  int ntpts; //number of tevol points
  int npoly;
  double field;
  double o1;
  
  // analysis specificaitons
  Matrix compspec;

  //inference options
  bool inferdw;
  bool inferk;

  //other options
  bool ffton;
  
  // ard flags
  bool doard;
 
};
