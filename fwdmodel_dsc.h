/*  fwdmodel_asl_grase.h - Implements the GRASE model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class DSCFwdModel : public FwdModel {
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
  { return 2 + (infermtt?1:0) + (inferlambda?1:0) + (inferdelay?1:0) + (inferart?2:0) + (inferret?1:0)+ (usecbv?1:0) + (dispoption?2:0); } 

  virtual ~DSCFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  // Constructor
  DSCFwdModel(ArgsType& args);

protected: 

  ColumnVector aifshift( const ColumnVector& aif, const float delta, const float hdelt ) const;
  void createconvmtx( LowerTriangularMatrix& A, const ColumnVector aifnew ) const;
  
// Constants

  // Lookup the starting indices of the parameters
  int cbf_index() const {return 1;} 

  int gmu_index() const {  return 1 + (infermtt?1:0);  }

  int lambda_index() const { return 1 + (infermtt?1:0) + (inferlambda?1:0); }

  int delta_index() const { return 1 + (infermtt?1:0) + (inferlambda?1:0) + (inferdelay?1:0); }
 
  int sig0_index() const { return 2 + (infermtt?1:0) + (inferlambda?1:0) + (inferdelay?1:0); }

  int art_index() const { return sig0_index() + (inferart?1:0);}

  int ret_index() const { return art_index() + (inferart?1:0) + (inferret?1:0); } //NB two arterial parameters

  int cbv_index() const { return ret_index() + (usecbv?1:0); }

  int disp_index() const { return cbv_index() + (dispoption?1:0); }

  //for ARD
  vector<int> ard_index;

    // scan parameters
  double te;
  double r2;
  double delt;
  ColumnVector artsig;
  ColumnVector s;

  bool aifconc;

  bool infermtt;
  bool usecbv;
  bool inferlambda;
  bool inferdelay;
  bool inferart;
  bool artoption;
  bool dispoption;
  bool inferret;
  bool doard;

  bool imageprior;

  string convmtx;

};
