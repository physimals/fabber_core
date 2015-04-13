/*  fwdmodel_asl_multiphase.h - 

    Michael Chappell, QuBIc (IBME) & FMRIB Image Analysis Group

    Copyright (C) 2013 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class MultiPhaseASLFwdModel : public FwdModel {
public: 
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;
  virtual void Initialise(MVNDist& posterior) const;
  static void ModelUsage();
  virtual string ModelVersion() const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const 
  { return 3 + (incvel?1:0) ; } 

  virtual ~MultiPhaseASLFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  // Constructor
  MultiPhaseASLFwdModel(ArgsType& args);


protected: // Constants
  int repeats;

  // modulation function
  string modfn;

  // inference options
  bool incvel;
  bool infervel;

  // fermi function variables
  double alpha;
  double beta;

  //modulation matrix
  double mod_fn( const double inphase, const double v) const;
  double interp(const ColumnVector& x, const ColumnVector& y, const double xi) const;
  Matrix mod_mat;
  ColumnVector mod_phase;
  ColumnVector mod_v;
  double vmax;
  double vmin;
  int nvelpts;

};
