/*  fwdmodel_asl_grase.h - Implements the GRASE model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#ifndef __FABBER_ASL_GRASE_FWDMODEL_H
#define __FABBER_ASL_GRASE_FWDMODEL_H 1

#include "fwdmodel.h"
#include "inference.h"
#include <string>

using namespace std;

class GraseFwdModel : public FwdModel {
public: 
  static FwdModel* NewInstance();

  // Virtual function overrides
  virtual void Initialize(ArgsType& args);
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;
  virtual vector<string> GetUsage() const;
  virtual string ModelVersion() const;
                  
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const 
  { return 2 - (singleti?1:0) + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0);

    //return 2 - (singleti?1:0) + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (inferinveff?1:0) + (infertrailing?1:0) + (infertaub?1:0); 
  } 

  virtual ~GraseFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;


protected: // Constants

  // Lookup the starting indices of the parameters
  int tiss_index() const {return 1;} //main tissue parameters: ftiss and delttiss alway come first

  int tau_index() const {  return 2 + (infertau?1:0);  }

  int art_index() const {  return 2 + (infertau?1:0) + (inferart?1:0); }

  int t1_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?1:0); }
  
  //int inveff_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) +(inferinveff?1:0); }

  //int trailing_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (infertrailing?1:0); }

  //int taub_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (inferinveff?1:0) + (infertrailing?1:0) + (infertaub?1:0);}

  int taub_index() const { return 2 + (infertau?1:0) + (inferart?2:0) + (infert1?2:0) + (infertaub?1:0);}

  // index for the parameter to expereicne ARD (this is the arterial perfusion flow)
  int ard_index() const { return 2 + (infertau?1:0) + (inferart?1:0); }
  
  // scan parameters
  double seqtau; //bolus length as set by the sequence
  double setdelt; //BAT for prior (tissue compartment)
  double deltprec; //precision for BAT
  int repeats;
  double t1;
  double t1b;
  double lambda;
  double pretisat;
  double slicedt;
  bool casl;
  bool grase; //to indicate data was collected with GRASE-ASL
  bool calib; //indicates calibrated data for T_1app calculation
  bool singleti; //specifies that only tissue perfusion should be inferred
  bool infertau;
  bool infertaub;
  bool inferart;
  bool infert1;
  //bool inferinveff;
  //bool infertrailing;
  bool doard;
  ColumnVector tis;
  Real timax;

 private:
  /** Auto-register with forward model factory. */
  static FactoryRegistration<FwdModelFactory, GraseFwdModel> registration;
  

};

#endif  // __FABBER_ASL_GRASE_FWDMODEL_H
