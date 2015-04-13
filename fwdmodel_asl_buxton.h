/*  fwdmodel_asl_buxton.h - Implements the Buxton kinetic curve model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class BuxtonFwdModel : public FwdModel {
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
  { return 2 + (infertau?1:0) + (infert1?2:0) + (twobol?2:0); } 

  virtual ~BuxtonFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

 virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;

  // Constructor
  BuxtonFwdModel(ArgsType& args);


protected: // Constants

  // Lookup the starting indices of the parameters
  int tiss_index() const {return 1;} //main tissue parameters: ftiss and delttiss always come first
  int tau_index() const { 
    if (infertau) { return 3; }
    else {return 0; } //zero mean parameter not set
  }
  int t1_index() const {
    if (infert1) {
      return 2 + (infertau?1:0) + 1;
    }
    else { return 0; }
  }
  int tiss2_index() const 
  {
    if (twobol) {
      return 2 + (infertau?1:0) + (infert1?2:0) +1;
    }
    else { return 0; }
  }

  // index for parameter subject to ARD (on perfusion of second bolus)
  int ard_index() const { return 2 + (infertau?1:0) + (infert1?2:0) + 1; }
  
  // scan parameters
  double seqtau; //bolus length as set by the sequence
  int repeats;
  double t1;
  double t1b;
  double lambda;
  bool infertau;
  bool infert1;
  bool twobol;
  bool doard;
  ColumnVector tis;
  Real timax;


};
