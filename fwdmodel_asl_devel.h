/*  fwdmodel_asl_devel.h - Development resting state ASL model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2009 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class DevelFwdModel : public FwdModel {
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
  { return (infertiss?2:0) - (singleti?1:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?(2+(infertau?1:0)+(infert1?1:0)+(usepve?2:0)):0)+2;  

    //return 2 - (singleti?1:0) + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (inferinveff?1:0) + (infertrailing?1:0) + (infertaub?1:0); 
  } 

  virtual ~DevelFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard);
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;

  // Constructor
  DevelFwdModel(ArgsType& args);


protected: // Constants

  // Lookup the starting indices of the parameters
  int tiss_index() const {return (infertiss?1:0);} //main tissue parameters: ftiss and delttiss alway come first

  int tau_index() const {  return (infertiss?2:0) + (infertiss?(infertau?1:0):0);  }

  int art_index() const {  return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?1:0); }

  int t1_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?1:0); }
  
  //int inveff_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) +(inferinveff?1:0); }

  //int trailing_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (infertrailing?1:0); }

  //int taub_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (inferinveff?1:0) + (infertrailing?1:0) + (infertaub?1:0);}

  int taub_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0);}

  //int R_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0) + (inferart?1:0);}
  int wm_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?1:0); }

  int pv_index() const { return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?(2 + (infertau?1:0) + (infert1?1:0) ):0) + (usepve?1:0); }

  int disp_index() const {return (infertiss?2:0) + (infertiss?(infertau?1:0):0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0)  + (inferwm?(2 + (infertau?1:0) + (infert1?1:0) ):0) + (usepve?2:0) +1; }

  // vector indices for the parameters to expereicne ARD
  vector<int> ard_index;


  // scan parameters
  double seqtau; //bolus length as set by the sequence
  int repeats;
  double t1;
  double t1b;
  double t1wm;
  double lambda;
  double pretisat;
  bool grase; //to indicate data was collected with GRASE-ASL

  bool infertiss;
  bool singleti; //specifies that only tissue perfusion should be inferred
  bool infertau;
  bool infertaub;
  bool inferart;
  bool infert1;
  bool inferwm;
  bool usepve;
  //bool inferinveff;
  //bool infertrailing;

  // ard flags
  bool doard;
  bool tissard;
  bool artard;
  bool wmard;

  ColumnVector tis;
  Real timax;

  //disperison parameters (that need to be availibe for icgf)
  float icgf(float a, float x) const;

};
