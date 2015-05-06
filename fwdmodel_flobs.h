/*  fwdmodel_flobs.h - Does FLOBS

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */


#ifndef __FABBER_FWDMODEL_FLOBS_H
#define __FABBER_FWDMODEL_FLOBS_H 1

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class FlobsFwdModel : public FwdModel {
public:
  static FwdModel* NewInstance();
public: 
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const
  {};
              
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const
  { assert(basis.Ncols()>0); 
    return basis.Ncols() + nuisanceBasis.Ncols(); }
  
  virtual string ModelVersion() const;

  virtual ~FlobsFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  virtual void Initialize(ArgsType& args);

protected: // Constants

  Matrix basis; // Ntr by Nbasis

  // Should also have a nuisance basis -- e.g. for offset
  Matrix nuisanceBasis; // Ntr by Nnuisance
};

// flobs7 = polar coordinate parameterization
// b1 cos(b2) X1 + b1 sin(b2) X2
// prior on b2 should probably be modified slightly (unless it's very small).
class Flobs7FwdModel : public FlobsFwdModel {
  public:
    static FwdModel* NewInstance();
    virtual vector<string> GetUsage() const;
    virtual void Evaluate(const ColumnVector& params, 
                          ColumnVector& result) const;
    virtual ~Flobs7FwdModel() { return; }

private:
  /** Auto-register with forward model factory. */
  static FactoryRegistration<FwdModelFactory, Flobs7FwdModel> registration;
};


// flobs5 = b1X1 + b1b2X2 parameterization
class Flobs5FwdModel : public FlobsFwdModel {
  public:
    static FwdModel* NewInstance();
    virtual vector<string> GetUsage() const;
    virtual void Evaluate(const ColumnVector& params,  
                          ColumnVector& result) const;
    virtual ~Flobs5FwdModel() { return; }

private:
  /** Auto-register with forward model factory. */
  static FactoryRegistration<FwdModelFactory, Flobs5FwdModel> registration;
};

#endif /* __FABBER_FWDMODEL_FLOBS_H */

