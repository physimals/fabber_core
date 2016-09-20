/*  fwdmodel.h - The base class for generic forward models

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

/* fwdmodel.h
 * Class declaration for generic forward models and related classes.
 * Written by Adrian Groves, 2007
 * FMRIB Centre, University of Oxford
 *
 * Last modified: $Date: 2013/04/29 12:38:19 $ $Author: chappell $ $Revision: 1.20 $
 */

#ifndef __FABBER_FWDMODEL_H
#define __FABBER_FWDMODEL_H

#include "assert.h"
#include <string>
#include <vector>

#include "newmatap.h"

#include "dist_mvn.h"
#include "utils.h"

using namespace std;
using namespace NEWMAT;

class FwdModel
{
public:

	/**
	 * Static member function, to pick a forward model from a name
	 */
	static FwdModel* NewFromName(const string& name);

	/**
	 * Get usage information for a named model
	 */
	static void UsageFromName(const string& name, std::ostream &stream);

	/**
	 * Initialize a new instance using configuration from the given
	 * arguments.
	 * @param args Configuration parameters.
	 */
	virtual void Initialize(FabberRunData& args) = 0;

	/**
	 * Evaluate the forward model
	 */
	virtual void Evaluate(const ColumnVector& params, ColumnVector& result) const = 0;

	/**
	 * Evaluate the gradient, the int return is to indicate whether a valid gradient is returned by the model
	 */
	virtual int Gradient(const ColumnVector& params, Matrix& grad) const;

	/**
	 * Return model usage information.
	 * @return vector of strings, one per line of information.
	 */
	virtual void Usage(std::ostream &stream) const;

	/**
	 * See fwdmodel.cc for an example of how to implement this.
	 *
	 * @return a CVS version info string
	 */
	virtual string ModelVersion() const;

	/**
	 * How many parameters in the model
	 */
	virtual int NumParams() const = 0;

	/**
	 * How many outputs for the model
	 *
	 * The default implementation calls Evaluate() on some fake data
	 * and sees how many outputs are produced.
	 */
	virtual int NumOutputs() const;

	/**
	 * Load up some sensible suggestions for initial prior & posterior values
	 *
	 * @param prior Prior distribution for parameters. Should be passed in as
	 *        an MVN of the correct size for the number of parameters. Model
	 *        may set mean/variances to suggested prior, or just give a trivial
	 *        default
	 * @param posterior As above for posterior distribution
	 */
	virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const = 0;

	/**
	 * Voxelwise initialization of the posterior, i.e. a parameter initialisation
	 *
	 * Called for each voxel, rather than HardcodedInitialDists which is called
	 * at the very start. Does not need to do anything unless the model
	 * wants per-voxel default posterior.
	 *
	 * @param posterior Posterior distribution for model parameters. Should be passed in as
	 *        an MVN of the correct size for the number of parameters. Model
	 *        may set mean/variances to suggested prior, or just do nothing to
	 *        accept general default
	 */
	virtual void InitParams(MVNDist& posterior) const
	{
	}

	/**
	 * Name each of the parameters
	 *
	 * See fwdmodel_linear.h for a generic implementation
	 */
	virtual void NameParams(vector<string>& names) const = 0;

	/**
	 * Describe what a given parameter vector means (to LOG)
	 *
	 * Default implementation uses NameParams to give reasonably meaningful output
	 */
	virtual void DumpParameters(const ColumnVector& params, const string& indent = "") const;

	/**
	 * An ARD update step can be specified in the model
	 */
	virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const
	{
	}
	;

	/**
	 * Setup function for the ARD process
	 *
	 * Forces the prior on the parameter that is subject to ARD to be correct -
	 * really a worst case scenario if people are loading in their own priors
	 */
	virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const
	{
	}
	;

	/**
	 * Indicies of parameters to which ARD should be applied
	 */
	vector<int> ardindices;

	/**
	 * For models that need the data values in the voxel to calculate
	 *
	 * Called for each voxel so the model knows what the current data is
	 */
	virtual void pass_in_data(const ColumnVector& voxdata)
	{
		data = voxdata;
	}

	/**
	 * For models that need the data and supplementary data values in the voxel to calculate
	 *
	 * Called for each voxel so the model knows what the current data is
	 */
	virtual void pass_in_data(const ColumnVector& voxdata, const ColumnVector& voxsuppdata)
	{
		data = voxdata;
		suppdata = voxsuppdata;
	}

	/**
	 * For models that need to know the voxel co-ordinates of the data
	 *
	 * Called for each voxel so the model knows what the current coords are
	 */
	virtual void pass_in_coords(const ColumnVector& coords);

	virtual ~FwdModel()
	{
	}

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

/** 
 * \ref SingletonFactory that returns pointers to \ref FwdModel.
 */
typedef SingletonFactory<FwdModel> FwdModelFactory;

#endif /* __FABBER_FWDMODEL_H */

