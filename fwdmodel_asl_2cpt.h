/*  fwdmodel_asl_2cpt.h - Implements the 2cpt model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2009 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class TwoCptFwdModel : public FwdModel {
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
  { return 2 - (singleti?1:0) + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0) + (inferPS?1:0);

    //return 2 - (singleti?1:0) + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (inferinveff?1:0) + (infertrailing?1:0) + (infertaub?1:0); 
  } 

  virtual ~TwoCptFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;

  // Constructor
  TwoCptFwdModel(ArgsType& args);


protected: // Constants

  // Lookup the starting indices of the parameters
  int tiss_index() const {return 1;} //main tissue parameters: ftiss and delttiss alway come first

  int tau_index() const {  return 2 + (infertau?1:0);  }

  int art_index() const {  return 2 + (infertau?1:0) + (inferart?1:0); }

  int t1_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?1:0); }
  
  int taub_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0);}

  int PS_index() const { return taub_index() + 1; }

  // index for the parameter to expereicne ARD (this is the arterial perfusion flow)
  int ard_index() const { return 2 + (infertau?1:0) + (inferart?1:0); }
  
  // scan parameters
  double seqtau; //bolus length as set by the sequence
  int repeats;
  double t1;
  double t1b;
  double pretisat;
  bool grase; //to indicate data was collected with GRASE-ASL
  bool calib; //indicates calibrated data for T_1app calculation
  bool slow; //true to use 'slow solution'
  bool singleti; //specifies that only tissue perfusion should be inferred
  bool infertau;
  bool infertaub;
  bool inferart;
  bool infert1;
  bool inferPS;
  //bool inferinveff;
  //bool infertrailing;
  bool doard;
  ColumnVector tis;
  Real timax;


};
