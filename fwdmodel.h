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

#include "utils.h"
#include "dist_mvn.h"

#include "newmat.h"

#include <string>
#include <vector>

class FwdModel
{
public:

	static void LoadFromDynamicLibrary(std::string filename);

	/**
	 * Static member function to return the names of all known
	 * models
	 */
	static std::vector<std::string> GetKnown();

	/**
	 * Static member function, to pick a forward model from a name
	 */
	static FwdModel* NewFromName(const std::string& name);

	/**
	 * Get usage information for a named model
	 */
	static void UsageFromName(const std::string& name, std::ostream &stream);

	/**
	 * Get option descriptions for this model. The default returns
	 * nothing to enable compatibility with older model code
	 */
	virtual void GetOptions(std::vector<OptionSpec> &opts) const {}
	
	/**
	 * @return human-readable description of the model.
	 */
	virtual std::string GetDescription() const {return "";}

	/**
	 * Get the model version. There is no fixed format for this,
	 * and it has no meaning other than by comparison with different
	 * versions of the same model. 
	 * 
	 * See fwdmodel.cc for an example 
	 * of how to implement this to return a CVS file version.
	 *
	 * @return a string indicating the model version. 
	 */
	virtual std::string ModelVersion() const;

	/**
	 * Initialize a new instance using configuration from the given
	 * arguments.
	 * 
	 * @param args Configuration parameters.
	 */
	virtual void Initialize(FabberRunData& args) = 0;

	/**
	 * How many parameters in the model? 
	 * 
	 * Initialize must be called before this method
	 * 
	 * @return number of parameters, i.e. size of vector to be passed
	 * to Evaluate function
	 */
	virtual int NumParams() const = 0;

	/**
	 * Name each of the parameters
	 * 
	 * Initialize must be called before this method
	 *
	 * See fwdmodel_linear.h for a generic implementation
	 */
	virtual void NameParams(std::vector<std::string>& names) const = 0;
	
	/**
	 * Load up some sensible suggestions for initial prior & posterior values
	 * 
	 * Initialize must be called before this method
	 *
	 * @param prior Prior distribution for parameters. Should be passed in as
	 *        an MVN of the correct size for the number of parameters. Model
	 *        may set mean/variances to suggested prior, or just give a trivial
	 *        default
	 *        
	 * @param posterior As above for posterior distribution
	 */
	virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const = 0;

	/**
	 * For models that need the data values in the voxel to calculate
	 *
	 * Called for each voxel so the model knows what the current data is
	 *
	 * @param voxdata Vector containing current voxel data. Evaluate will return the same
	 *                number of values that are in this vector
	 */
	virtual void pass_in_data(const NEWMAT::ColumnVector& voxdata)
	{
		data = voxdata;
	}

	/**
	 * For models that need the data and supplementary data values in the voxel to calculate
	 *
	 * Called for each voxel so the model knows what the current data is
	 * 
	 * @param voxdata Vector containing current voxel data. Evaluate will return the same
	 *                number of values that are in this vector
	 * @param voxsuppdata Supplementary data if provided. Must be the same length as voxdata.
	 */
	virtual void pass_in_data(const NEWMAT::ColumnVector& voxdata, const NEWMAT::ColumnVector& voxsuppdata)
	{
		data = voxdata;
		suppdata = voxsuppdata;
	}

	/**
	 * For models that need to know the voxel co-ordinates of the data
	 *
	 * Called for each voxel so the model knows what the current coords are
	 * 
	 * @param coords Vector of length 3 containing x, y, z coords
	 */
	virtual void pass_in_coords(const NEWMAT::ColumnVector& coords);

	/**
	 * Voxelwise initialization of the posterior, i.e. a parameter initialisation
	 *
	 * Called for each voxel, rather than HardcodedInitialDists which is called
	 * at the very start. Does not need to do anything unless the model
	 * wants per-voxel default posterior.
	 * 
	 * Initialize must be called before this method
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
	 * Evaluate the forward model
	 * 
	 * Initialize must be called before this method
	 * 
	 * @param params Model parameter values. Must contain the correct number of parameters
	 *  			 as specified by NumParams
	 * @param result Will be populated with the model prediction for these parameters.
	 *               The length of this vector will be set to the same as the number of
	 *               data points passed in via pass_in_data
	 */
	virtual void Evaluate(const NEWMAT::ColumnVector& params, NEWMAT::ColumnVector& result) const = 0;

	/**
	 * Evaluate the gradient
	 * 
	 * @param params Model parameter values. Must contain the correct number of parameters
	 *  			 as specified by NumParams
	 * @param grad If returning true, gradient of model as a NumParams x NumParams matrix
	 * 
	 * @return false if no valid gradient is returned by the model, true if it is
	 */
	virtual bool Gradient(const NEWMAT::ColumnVector& params, NEWMAT::Matrix& grad) const;

	/**
	 * An ARD update step can be specified in the model
	 */
	virtual void UpdateARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const
	{
	}
	
	/**
	 * Setup function for the ARD process
	 *
	 * Forces the prior on the parameter that is subject to ARD to be correct -
	 * really a worst case scenario if people are loading in their own priors
	 */
	virtual void SetupARD(const MVNDist& posterior, MVNDist& prior, double& Fard) const
	{
	}

	/**
	 * Indicies of parameters to which ARD should be applied
	 */
	vector<int> ardindices;

	virtual ~FwdModel()
	{
	}

#ifdef DEPRECATED
	/**
	 * Describe what a given parameter vector means (to LOG)
	 *
	 * Default implementation uses NameParams to give reasonably meaningful output
	 */
	virtual void DumpParameters(const NEWMAT::ColumnVector& params, const std::string& indent = "") const;

	/**
	 * Return model usage information.
	 *
	 * Deprecated in favour of GetOptions.
	 *
	 * @return vector of strings, one per line of information.
	 */
	virtual void Usage(std::ostream &stream) const;
#endif
	
protected:
	// Your derived classes should have storage for all constants that are
	// implicitly part of g() -- e.g. pulse sequence parameters, any parameters
	// that are assumed to take known values, and basis functions.  Given these
	// constants, NumParams() should have a fixed value.
	
	// storage for voxel co-ordinates
	int coord_x;
	int coord_y;
	int coord_z;

	//storage for data
	NEWMAT::ColumnVector data;
	NEWMAT::ColumnVector suppdata;
};

/** 
 * \ref SingletonFactory that returns pointers to \ref FwdModel.
 */
typedef SingletonFactory<FwdModel> FwdModelFactory;

/**
 * Function pointer for a function which creates new instances of FwdModel.
 *
 * Used for the dynamic loading of models
 */
typedef FwdModel* (*NewInstanceFptr)(void);

#endif /* __FABBER_FWDMODEL_H */

