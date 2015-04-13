/*  fwdmodel_linear.h - Linear forward model and related classes

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "fwdmodel.h"

class LinearFwdModel : public FwdModel {
 public:
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
                      ColumnVector& result) const;
  virtual int NumParams() const { return centre.Nrows(); }
  virtual void DumpParameters(const ColumnVector& vec,
                              const string& indent = "") const;                            
  virtual void NameParams(vector<string>& names) const;

  ReturnMatrix Jacobian() const { return jacobian; }
  ReturnMatrix Centre() const { return centre; }
  ReturnMatrix Offset() const { return offset; }

  LinearFwdModel(const Matrix& jac, 
		 const ColumnVector& ctr, 
		 const ColumnVector& off) 
    : jacobian(jac), centre(ctr), offset(off) 
    { assert(jac.Nrows() == ctr.Ncols()); assert(jac.Ncols() == off.Ncols()); }
    
  // Upgrading to a full externally-accessible model type
  LinearFwdModel(ArgsType& args);
  virtual string ModelVersion() const;
  static void ModelUsage();
  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

 protected:
  LinearFwdModel() { return; } // Leave uninitialized; derived classes only

  Matrix jacobian;     // J (tranposed?)
  ColumnVector centre; // m
  ColumnVector offset; // g(m)
    // The amount to effectively subtract from Y is g(m)-J*m
};

class LinearizedFwdModel : public LinearFwdModel {
public:
  // Virtual function overrides
  virtual void DumpParameters(const ColumnVector& vec,
                              const string& indent = "") const;
  virtual void NameParams(vector<string>& names) const
    { assert(fcn); fcn->NameParams(names); }
  using FwdModel::ModelVersion; // Tell the compiler we want both ours and the base version
  string ModelVersion() { assert(fcn != NULL); return fcn->ModelVersion(); }

  // Constructor (leaves centre, offset and jacobian empty)
  LinearizedFwdModel(const FwdModel* model) : fcn(model) { return; }
  
  // Copy constructor (needed for using vector<LinearizedFwdModel>)
  // NOTE: This is a reference, not a pointer... and it *copies* the
  // given LinearizedFwdModel, rather than using it as its nonlinear model!
  LinearizedFwdModel(const LinearizedFwdModel& from) 
    : LinearFwdModel(from), fcn(from.fcn) { return; }

  void ReCentre(const ColumnVector& about);
  // centre=about; offset=fcn(about); 
  // jacobian = numerical differentiation about centre

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const
    { assert(fcn); fcn->HardcodedInitialDists(prior, posterior); }

  
private:
  const FwdModel* fcn;  
};

