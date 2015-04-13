/*  fwdmodel_fasl.h - multi-TI functional ASL model

    Michael Chappell, IBME & FMRIB Image Analysis Group

    Copyright (C) 20010 University of Oxford  */

/*  CCOPYRIGHT */


#ifndef __FABBER_FWDMODEL_FASL_H
#define __FABBER_FWDMODEL_FASL_H 1

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class FASLFwdModel : public FwdModel {
public: 
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;
                  
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const
  { return D0index()+Dbasis.Ncols() + ((stdevInvEff>0)?1:0) + ((stdevT1b>0)?1:0) + ((stdevT1>0)?1:0) +1;}
      // { return NnStart(echoTime.Ncols()+1) + (stdevT1b>0?1:0) 
    //   + (stdevInvEff>0?1:0) + (stdevDt>0?1:0) - 1; }

  virtual string ModelVersion() const;

  virtual ~FASLFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  // Constructor
  FASLFwdModel(ArgsType& args);
  // Usage info
  static void ModelUsage();


protected: // Constants
  
  float kctissue_nodisp(const float ti, const float delttiss, const float tau, const float T_1b, const float T_1app) const;
  
  int Q0index() const { return 1; }
  int M0index() const { return Q0index()+Qbasis.Ncols()+1; }
  int D0index() const {return M0index()+Mbasis.Ncols()+1; }
  //  int R0index() const { return M0index()+Mbasis.Ncols()+1; }
  //int NnStart(int te) const 
  //  { assert(te > 0 && te <= echoTime.Nrows()+1); // go one past, for getting end point
  //    return R0index()+Rbasis.Ncols()+(te-1)*Nbasis.Ncols() + 1; }
  int InvEffIndex() const { assert(stdevInvEff>0); return D0index()+Dbasis.Ncols() +1 ; }
  int T1bIndex() const { assert(stdevT1b>0); 
    return D0index()+Dbasis.Ncols()+((stdevInvEff>0)?1:0)+1; }
 int T1Index() const { assert(stdevT1>0); 
    return D0index()+Dbasis.Ncols()+((stdevInvEff>0)?1:0)+((stdevT1b>0)?1:0)+1; }
 int AIndex() const { 
   return D0index()+Dbasis.Ncols()+((stdevInvEff>0)?1:0)+((stdevT1b>0)?1:0)+((stdevT1>0)?1:0)+1; }
  //int dtIndex() const { assert(stdevDt>0);
  //  return NnStart(echoTime.Ncols()+1) + (stdevT1b>0?1:0) + (stdevInvEff>0?1:0); }

  // Slow submatrixing... but it works.
  // In fact, it's plenty fast enough.. fwd model is a small part calculation time
  ReturnMatrix QnOf(const ColumnVector& params) const
    { return params.Rows(Q0index()+1,Q0index()+Qbasis.Ncols()).Evaluate(); }
 ReturnMatrix DnOf(const ColumnVector& params) const
    { return params.Rows(D0index()+1,D0index()+Dbasis.Ncols()).Evaluate(); }
  ReturnMatrix MnOf(const ColumnVector& params) const
    { return params.Rows(M0index()+1,M0index()+Mbasis.Ncols()).Evaluate(); }
//  ReturnMatrix RnOf(const ColumnVector& params) const
//    { return params.Rows(R0index()+1,R0index()+Rbasis.Ncols()).Evaluate(); }  
//  ReturnMatrix NnOf(int te, const ColumnVector& params) const
//    { return params.Rows(NnStart(te),NnStart(te+1)-1).Evaluate(); }

  // scan parameters
  ColumnVector rho;
//ColumnVector echoTime;
ColumnVector tivec;
float tau;
  
 bool dataisdiff;

  // assumed known constants
 double fixedInvEff, fixedT1b, fixedDt, fixedT1;
 double stdevT1b, stdevInvEff, stdevDt, stdevT1; // >0 means treat these as parameters!
  
  // basis for each piece
  Matrix Qbasis; // basis in columns?
Matrix Dbasis;
  Matrix Mbasis;
//  Matrix Rbasis;
//  Matrix Nbasis;
};

#endif /* __FABBER_FWDMODEL_QUIPSS2_H */

