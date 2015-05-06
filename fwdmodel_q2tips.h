/*  fwdmodel_q2tips.h - Implements the Q2TIPS ASL model

    Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & QuBIc (IBME)

    Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_quipss2.h"

// The Q2TIPS model is almost identical to QUIPSS II model.
// Only the Evaluate() function needs to change (and only slightly)

class Q2tipsFwdModel : public Quipss2FwdModel {

public: 
  static FwdModel* NewInstance();

  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;

  virtual string ModelVersion() const;

  virtual ~Q2tipsFwdModel() { return; }

  virtual void Initialize(ArgsType& args);

  virtual vector<string> GetUsage() const;

private:
  /** Auto-register with forward model factory. */
  static FactoryRegistration<FwdModelFactory, Q2tipsFwdModel> registration;

};
