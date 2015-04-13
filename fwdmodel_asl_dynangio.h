/*  fwdmodel_asl_quasar.h - Resting state ASL model for QUASAR acquisitions

    Michael Chappell, IBME & FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>

#include "asl_models.h"

using namespace std;
using namespace OXASL;

class DynAngioFwdModel : public FwdModel {
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
  { return 7;  
  } 

  virtual ~DynAngioFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard);
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;

  // Constructor
  DynAngioFwdModel(ArgsType& args);


protected: // Constants

  // Lookup the starting indices of the parameters

  int disp_index() const {return 4 +1 ; }


  // vector indices for the parameters to expereicne ARD
  vector<int> ard_index;


  // scan parameters
  double seqtau; //bolus length as set by the sequence;
  double t1b;
  float dti; //TI interval
  float FA; //flip angle

  string disptype;

  ColumnVector tis;
  Real timax;

  bool taufix;

};
