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
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const;
                  
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const
  { assert(basis.Ncols()>0); 
    return basis.Ncols() + nuisanceBasis.Ncols(); }
  
  virtual string ModelVersion() const;

  virtual ~FlobsFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  // Constructor
  FlobsFwdModel(ArgsType& args, bool sepScale) ;
  // Usage info
  static void ModelUsage();

protected: // Constants
 
  //  const bool useSeparateScale;
  const bool usePolarCoordinates;

  Matrix basis; // Ntr by Nbasis

  // Should also have a nuisance basis -- e.g. for offset
  Matrix nuisanceBasis; // Ntr by Nnuisance
};

#endif /* __FABBER_FWDMODEL_FLOBS_H */

