/*  fwdmodel_asl_rest.cc - Resting state ASL model

    Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

    Copyright (C) 2011-15 University of Oxford  */

/*  CCOPYRIGHT */

#ifndef __FABBER_ASL_REST_FWDMODEL_H
#define __FABBER_ASL_REST_FWDMODEL_H 1

#include "fwdmodel.h"
#include "inference.h"
#include <string>

#include "asl_models.h"

using namespace std;
using namespace OXASL;

class ASLFwdModel : public FwdModel {
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
  { return (inctiss?(incbat?2:1):0) + (incwm?(incbat?2:1):0) + (incart?(incbat?2:1):0) + (inctau?(septau?ncomps:1):0) + (inct1?(1+(inctiss?1:0)+(incwm?1:0)):0) + (incpve?2:0) + (incpc?( (inctiss?1:0) + (incwm?1:0)):0) + (incdisp?(sepdisp?( (inctiss?tiss_model->NumDisp():0) + (incwm?tiss_model->NumDisp():0) + (incart?art_model->NumDisp():0) ):tiss_model->NumDisp()):0) + (incexch?tiss_model->NumResid():0) + (incfacorr?1:0) + (incstattiss?1:0) ;
  } 

  virtual ~ASLFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;
  virtual void InitParams(MVNDist& posterior) const;

  virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;
  vector<int> ardindices;


protected: // Constants

  // Parameter Map:
  // Flow (GM, WM, art);
  // BAT (GM, WM, art);
  // tau (shared or GM,WM,art);
  // t1 (GM, WM, blood);
  // PV (GM, WM);
  // taupc (GM,WM);
  // disp (shared or GM,WM,art);

  // Lookup the starting indices of the parameters
  int flow_index() const {return 1;}

  int bat_index() const {return flow_index() + ncomps;}

  int tau_index() const {  return bat_index() + (incbat?ncomps:0);  }

  int t1_index() const { return tau_index() + (inctau?(septau?ncomps:1):0); } // special case since tau may not be separated
  
  int pv_index() const { return t1_index() + (inct1?( (inctiss?1:0) + (incwm?1:0) + 1 ):0); } //special case since we always have T1 of blood

  int taupc_index() const { return pv_index() + (incpve?2:0); }

  int disp_index() const { return taupc_index() + (incpc?( (inctiss?1:0) + (incwm?1:0) ):0); } // NB in the absence of any tissue compoennts (when this parameter is meaningless) this reduces to taupc_index and thus still makes sense in sequence

  int resid_index() const { return disp_index() +  (incdisp?(sepdisp?( (inctiss?tiss_model->NumDisp():0) + (incwm?tiss_model->NumDisp():0) + (incart?art_model->NumDisp():0) ):tiss_model->NumDisp()):0); }

  int facorr_index() const { return resid_index() + (inctiss?resid_model->NumResid():0); } //flip angle correction for LL - note that we assume that residue function parameters only come from the tissue (and WM shares the same ones if it has been added)

int stattiss_index() const {return facorr_index() + (incfacorr?1:0);} //including static tissue to the signal (e.g. for non-subtracted data)

  // index for the parameter to expereicne ARD (this is the arterial perfusion flow)
  int ard_index() const { return 1 + (inctiss?1:0) + (incwm?1:0); } // refers to aBV (art 'flow'), NB value still makes sense in the (illogical) condition of only arterial component and ARD
  
  // scan parameters
  double seqtau; //bolus length as set by the sequence
  double pretisat;
  double slicedt;
  bool casl;
  ColumnVector tis;
  ColumnVector taus; // bolus durations if there are more than one
  Real timax;
  int repeats;
  ColumnVector crush;
  Matrix crushdir;

  //data information
  bool raw; //data is 'raw' (non-subtracted)
  bool tagfirst; // raw data has tag (true) or control (false) first

  // analysis options/settings
  double setdelt; //BAT for prior (tissue compartment)
  double setdeltwm;
  double setdeltart;
  double deltprec; //precision for BAT
  bool calib; //indicates calibrated data for T_1app calculation
  bool doard;

  // T1 values
  double t1;
  double t1b;
  double t1wm;

  // Other model parameters
  double lambda;
  double lamgm;
  double lamwm;

  // Look-Locker
  bool looklocker;
  double FA; //the flip angle (requested fo the scanner)
  double dti; //inter TI spacing
  bool incfacorr;
  double dg; //flip angle correction term (for slice profile)

  // Hadamard time encoding
  bool hadamard; // indicates that we are modelling hadamard data
  int HadamardSize; // Size of the hadamard encoding
  int NumberOfSubBoli; // number of sub boli (might just different from HadamardSize)
  Matrix HadEncMatrix; //stores the Hadamard encoding matrix
  Matrix HadamardMatrix( const int size) const; // function to generate the hadamard matrix

  // inference/inclusion
  // -components
  bool inctiss; //tissue component (used for GM)
  bool infertiss;
  bool incart; //arterial (MV) component
  bool inferart;
  bool incwm; // white matter component
  bool inferwm;
  int ncomps; // number of components
  // -common things
  bool incbat; // bolus arrival time
  bool inferbat;
  bool incpc; //pre-capiliary (this can be GM and WM) - separated (different values for different components)
  bool inferpc;
  bool inctau; // bolus duration (art,GM,WM) - seperable
  bool infertau;
  bool septau; //separate tau values for different components
  bool multitau; //multiple tau for each TI
  bool inct1; // T1 values (blood,GM,WM)
  bool infert1;
  bool incdisp; // Dispersion (art,GM,WM) - separable
  bool inferdisp;
  bool sepdisp; //separate dispersion for different components
  // int ndisp; //number of dispersion parameters (per component)
  bool incexch; //restricted exchange (GM,WM) - separated
  bool inferexch;
  // -special
  bool incpve; //include partial volume estimates for PV correction (BUT do not use them unless pvcorr is set)
  bool pvcorr; //do PV correction
  // - static tissue
  bool incstattiss; // include a static tissue contribution to the signal (e.g. in case of non-subtracted data)
  bool inferstattiss; //

  //relative indices
  int wmidx;
  int artidx;

  // Models to use
  AIFModel* art_model;
  TissueModel* tiss_model;
  TissueModel* pc_model; //the model for pre-capilliary component
  ResidModel* resid_model;
  string disptype;
  string exchtype;

private:
  /** Auto-register with forward model factory. */
  static FactoryRegistration<FwdModelFactory, ASLFwdModel> registration;

};

#endif  //__FABBER_ASL_REST_FWDMODEL_H
