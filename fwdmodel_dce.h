/*  fwdmodel_dce.h - Implements the Dynamic Contrast Enhanced MRI model

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class DCEFwdModel : public FwdModel {
public: 
  // Virtual function overrides
  virtual void Evaluate(const ColumnVector& params, ColumnVector& result) const;
  //virtual void Evaluate(Volume Kev, Volume T10, Volume Ve, Volume oi3d) const;
  static void ModelUsage();
  virtual string ModelVersion() const;
                  
  virtual void DumpParameters(const ColumnVector& vec,
                                const string& indents = "") const;
                                
  virtual void NameParams(vector<string>& names) const;     
  virtual int NumParams() const 
  { return 5; } 

  virtual ~DCEFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

 virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;

  // Constructor
  DCEFwdModel(ArgsType& args);


protected: // Constants

  // Lookup the starting indices of the parameters

  // index for parameter subject to ARD (on perfusion of second bolus)
  int doard;
  int ard_index() const { return 1; }
  
  // scan parameters
  int tdim;  // number of time points

};
