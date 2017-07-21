/*  fwdmodel.h - The base class for generic forward models

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "dist_mvn.h"
#include "easylog.h"
#include "factories.h"
#include "transforms.h"

#include <newmat.h>

#include <string>
#include <vector>

/**
 * Information about a model parameter
 */
struct Parameter
{
	Parameter(unsigned int idx, std::string name="", DistParams prior=DistParams(0, 1), DistParams post=DistParams(0, 1), char prior_type='N', const Transform *transform=TRANSFORM_IDENTITY())
	  : idx(idx), name(name), prior(prior), post(post), prior_type(prior_type), transform(transform) {}
	
	unsigned int idx;
	std::string name;
	DistParams prior;
	DistParams post;
	char prior_type;
	const Transform  *transform;

	/** 
	 * Additional options, e.g. filename of data for image prior
	 *
	 * This should generally only by initialized by Fabber run options,
	 * not by the model itself 
	 */
	std::map<std::string, std::string> options; 
};

class FwdModel : public Loggable
{
public:

	/** Required in case subclasses manage resources */
    virtual ~FwdModel()
    {
    }

    /**
	 * Get option descriptions for this model. The default returns
	 * nothing to enable compatibility with older model code
	 */
    virtual void GetOptions(std::vector<OptionSpec> &opts) const {}

    /**
	 * @return human-readable description of the model.
	 */
    virtual std::string GetDescription() const;

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
	 * @param rundata Configuration parameters.
	 */
    virtual void Initialize(FabberRunData &rundata);

    /**
	 * Voxelwise initialization of the posterior in model space
	 *
	 * Subclasses can override this if they want to set the initial
	 * posterior on a per-voxel basis (e.g. an offset might be initially
	 * set to the value of the first data point).
	 * 
	 * Called by InitPostVox. This should set values in terms of model parameters
	 * and not worry about parameter transformations.
	 *
	 * @param posterior Initial posterior distribution for model parameters.
	 */
    virtual void InitParams(MVNDist &posterior) const
    {
    }
	
    /**
	 * Get parameter descriptions for this model. 
	 *
	 * The output is only valid after Initialize is called, as the
	 * parameters may well depend on options passed to the model
	 *
	 * For backwards compatibility, the default implementation 
	 * uses NameParams() and HardcodedInitialDists to supply this
	 * information, and assumes identity transforms
	 */
    void GetParameters(FabberRunData &rundata, std::vector<Parameter> &params);

    /**
	 * For models that need the data and supplementary data values in the voxel to calculate
	 *
	 * Called for each voxel so the model knows what the current data is
	 * 
	 * @param voxdata Vector containing current voxel data. Evaluate will return the same
	 *                number of values that are in this vector
	 * @param voxsuppdata Supplementary data if provided. Default is an empty vector. If not
	 *                    empty, must be the same length as voxdata.
	 */
    void PassData(const NEWMAT::ColumnVector &voxdata, const NEWMAT::ColumnVector &coords, const NEWMAT::ColumnVector &voxsuppdata=NEWMAT::ColumnVector());

    /**
	 * Initialization of the posterior.
	 *
	 * This is called for each voxel. The parameter defaults are used to set up 
	 * an initial posterior, then InitParams is called to allow the model to 
	 * do per-voxel initialization if required. Finally parameter transforms are
	 * applied so the resulting posterior contains appropriate values for Fabber's
	 * internal logic.
	 *
	 * Note that in doing the transformation, it is assumed that the initial posterior
	 * is diagonal.
	 * 
	 * Initialize must be called before this method
	 *
	 * @param posterior Posterior distribution for model parameters. Should be passed in as
	 *        an MVN of the correct size for the number of parameters. Model
	 *        may set mean/variances to suggested posterior, or just do nothing to
	 *        accept general default
	 */
    void GetInitialPosterior(MVNDist &posterior) const;

	/**
	 * Evaluate the forward model in Fabber internal space
	 * 
	 * This method calls the model-specific Evaluate method, but handles parameter transforms
	 * transparently.
	 *
	 * Initialize must be called before this method
	 * 
	 * @param params Model parameter values in Fabber internal space. 
	 * @param result Will be populated with the model prediction for these parameters.
	 *               The length of this vector will be set to the same as the number of
	 *               data points passed in via pass_in_data
	 */
	void EvaluateFabber(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;

    /**
	 * Evaluate the forward model in model space
	 * 
	 * Initialize must be called before this method
	 * 
	 * @param params Model parameter values. Must contain the correct number of parameters
	 *  			 as specified by NumParams
	 * @param result Will be populated with the model prediction for these parameters.
	 *               The length of this vector will be set to the same as the number of
	 *               data points passed in via pass_in_data
	 */
    virtual void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const = 0;

	/**
	 * Transform an MVN containing model values to Fabber internal values. 
	 *
	 * NB: Only the diagonal elements of the covariance are affected
	 */
	void ToFabber(MVNDist &mvn) const;

	/**
	 * Transform an MVN containing Fabber internal values to model values. 
	 *
	 * NB: Only the diagonal elements of the covariance are affected
	 */
	void ToModel(MVNDist &mvn) const;

	/**
	 * Load models from a dynamic library, adding them to the FwdModelFactory
	 */
    static void LoadFromDynamicLibrary(const std::string &filename, EasyLog *log = 0);

    /**
	 * Static member function to return the names of all known
	 * models
	 */
    static std::vector<std::string> GetKnown();

    /**
	 * Static member function, to pick a forward model from a name
	 */
    static FwdModel *NewFromName(const std::string &name);

    /**
	 * Get usage information for a named model
	 */
    static void UsageFromName(const std::string &name, std::ostream &stream);

#ifdef DEPRECATED
    /**
	 * How many parameters in the model? 
	 * 
	 * Initialize must be called before this method. DEPRECATED, replace with
	 * GetParameterDefaults()
	 * 
	 * @return number of parameters, i.e. size of vector to be passed
	 * to Evaluate function
	 */
    virtual int NumParams() const {return m_params.size();}

    /**
	 * Name each of the parameters
	 * 
	 * Initialize must be called before this method. DEPRECATED, replace with
	 * GetParameterDefaults()
	 */
    virtual void NameParams(std::vector<std::string> &names) const {}

    /**
	 * Load up some sensible suggestions for initial prior & posterior values
	 * 
	 * Initialize must be called before this method. DEPRECATED, replace with
	 * GetParameterDefaults()
	 *
	 * @param prior Prior distribution for parameters. Should be passed in as
	 *        an MVN of the correct size for the number of parameters. Model
	 *        may set mean/variances to suggested prior, or just give a trivial
	 *        default
	 *        
	 * @param posterior As above for posterior distribution
	 */
    virtual void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const {};

    /**
	 * An ARD update step can be specified in the model
	 */
    virtual void UpdateARD(const MVNDist &posterior, MVNDist &prior, double &Fard) const
    {
    }

    /**
	 * Setup function for the ARD process
	 *
	 * Forces the prior on the parameter that is subject to ARD to be correct -
	 * really a worst case scenario if people are loading in their own priors
	 */
    virtual void SetupARD(const MVNDist &posterior, MVNDist &prior, double &Fard) const
    {
    }

    /**
	 * Indicies of parameters to which ARD should be applied
	 */
    std::vector<int> ardindices;

    /**
	 * Describe what a given parameter vector means (to LOG)
	 *
	 * Default implementation uses NameParams to give reasonably meaningful output
	 */
    virtual void DumpParameters(const NEWMAT::ColumnVector &params, const std::string &indent = "") const;

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

    virtual void GetParameterDefaults(std::vector<Parameter> &params) const;

    // Your derived classes should have storage for all constants that are
    // implicitly part of g() -- e.g. pulse sequence parameters, any parameters
    // that are assumed to take known values, and basis functions.  Given these
    // constants, NumParams() should have a fixed value.

    // storage for current voxel data
    NEWMAT::ColumnVector coords;
    NEWMAT::ColumnVector data;
    NEWMAT::ColumnVector suppdata;

#ifdef DEPRECATED
	int coord_x;
	int coord_y;
	int coord_z;
#endif

private:
	std::vector<Parameter> m_params;
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
typedef FwdModel *(*NewInstanceFptr)(void);
