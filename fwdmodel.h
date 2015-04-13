/*  fwdmodel.h - The base class for generic forward models

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */


/* fwdmodel.h
 * Class declaration for generic forward models and related classes.
 * Written by Adrian Groves, 2007
 * FMRIB Centre, University of Oxford
 *
 * Last modified: $Date: 2013/04/29 12:38:19 $ $Author: chappell $ $Revision: 1.20 $
 */

//#pragma once // prevent multiple includes
#ifndef __FABBER_FWDMODEL_H
#define __FABBER_FWDMODEL_H 1

#include "assert.h"
#include "newmatap.h"
#include <string>
#include <vector>
#include "dist_mvn.h"
#include "easyoptions.h"

using namespace NEWMAT;
using namespace std;

/*
class FwdModelIdStruct {
public:
  virtual void DumpVector(const ColumnVector& vec, const string& indent = "") = 0;
  virtual ~FwdModelIdStruct() { return; }
};*/

class FwdModel {
public:
  // Virtual functions: common to all FwdModels
  
  virtual void Evaluate(const ColumnVector& params, 
			      ColumnVector& result) const = 0;
  // Evaluate the forward model

  virtual int Gradient(const ColumnVector& params, Matrix& grad) const;
  // evaluate the gradient, the int return is to indicate whether a valid gradient is returned by the model
                  
  virtual string ModelVersion() const; 
  // Return a CVS version info string
  // See fwdmodel.cc for an example of how to implement this.

  virtual int NumParams() const = 0;
  // How long should the parameter vector be?
  
  virtual int NumOutputs() const;
  // How long is output vector?  Default implementation uses Evaluate.

  // Various other useful functions:
  virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const = 0;
  // Load up some sensible suggestions for initial prior & posterior values

  virtual void Initialise(MVNDist& posterior) const {};
  // voxelwise initialization of the posterior
 
  virtual void NameParams(vector<string>& names) const = 0;
  // Name each of the parameters -- see fwdmodel_linear.h for a generic implementation
  
  virtual void DumpParameters(const ColumnVector& params, 
                              const string& indent="") const;
  // Describe what a given parameter vector means (to LOG)
  // Default implementation uses NameParams to give reasonably meaningful output 
  
  
  // Static member function, to pick a forward model from a name
  static FwdModel* NewFromName(const string& name, ArgsType& args);
  
  // Usage information for this model
  static void ModelUsageFromName(const string& name, ArgsType& args);

  // An ARD update step can be specified in the model
  virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const { return; };

  // Steup function for the ARD process (forces the prior on the parameter that is subject to ARD to be correct) - really a worst case scenario if people are loading in their own priors
  virtual void SetupARD( const MVNDist& posterior, MVNDist& prior, double& Fard ) const { return; };

  //vector of indicies of parameters to which ARD should be applied;
  vector<int> ardindices;
  
  // For models that need the data values in the voxel to calculate
  virtual void pass_in_data( const ColumnVector& voxdata ) { data=voxdata; return; };
  virtual void pass_in_data( const ColumnVector& voxdata, const ColumnVector& voxsuppdata )
  { data = voxdata; suppdata = voxsuppdata; return; };

  // For models that need to know the voxel co-ordinates of the data
  virtual void pass_in_coords( const ColumnVector& coords);

  virtual ~FwdModel() { return; };
  // Virtual destructor
  
  // Your derived classes should have storage for all constants that are
  // implicitly part of g() -- e.g. pulse sequence parameters, any parameters
  // that are assumed to take known values, and basis functions.  Given these
  // constants, NumParams() should have a fixed value.

 protected:
  // storage for voxel co-ordinates
  int coord_x;
  int coord_y;
  int coord_z;
  //storage for data
  ColumnVector data;
  ColumnVector suppdata;
};

#endif /* __FABBER_FWDMODEL_H */

