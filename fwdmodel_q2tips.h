/*  fwdmodel_q2tips.h - Implements the Q2TIPS ASL model

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_quipss2.h"

// The Q2TIPS model is almost identical to QUIPSS II model.
// Only the Evaluate() function needs to change (and only slightly)

class Q2tipsFwdModel : public Quipss2FwdModel {

public: 
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;

  virtual string ModelVersion() const;

  virtual ~Q2tipsFwdModel() { return; }

  // Constructor
  Q2tipsFwdModel(ArgsType& args) : Quipss2FwdModel(args) { }

};
