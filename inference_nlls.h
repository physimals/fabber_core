/* inference_nlls.h - Non-Linear Least Squares class declarations

   Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

   Copyright (C) 2007-2015 University of Oxford */

/*  CCOPYRIGHT  */

#include "inference.h"
#include "miscmaths/nonlin.h"
#include <boost/shared_ptr.hpp>
#include "miscmaths/bfmatrix.h"

class NLLSInferenceTechnique : public InferenceTechnique {
 public:
  static InferenceTechnique* NewInstance();
  NLLSInferenceTechnique() { return; }
  virtual void Initialize(FwdModel* fwd_model, ArgsType& args);
  virtual void DoCalculations(const DataSet& data);
  virtual ~NLLSInferenceTechnique();
 protected:
  const MVNDist* initialFwdPosterior;
  bool vbinit;
  bool lm;
};

class NLLSCF : public NonlinCF
{
 public:
 NLLSCF(const ColumnVector& pdata, const FwdModel* pm) 
   : y(pdata), model(pm), linear(pm) {}
  ~NLLSCF() { return; }
  virtual double cf(const ColumnVector& p) const;
  virtual ReturnMatrix grad(const ColumnVector& p) const;
  virtual boost::shared_ptr<BFMatrix> hess(const ColumnVector& p, boost::shared_ptr<BFMatrix> iptr) const;
 private:
  const ColumnVector y; //Values from data
  const FwdModel* model;
  mutable LinearizedFwdModel linear;
};
