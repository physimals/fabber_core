/*  fwdmodel_biexp.h - 

    Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class BiExpFwdModel : public FwdModel {
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
  { return 4; } 

  virtual ~BiExpFwdModel() { return; }

  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

  virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const;

  // Constructor
  BiExpFwdModel(ArgsType& args);


protected: // Constants

  
  // index for the parameter to expereicne ARD 
  int ard_index() const { return 3; }
  
  // scan parameters
 
  bool doard;
  

};
