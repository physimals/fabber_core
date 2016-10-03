/*  noisemodel.h - Class declaration for generic noise models

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#ifndef __FABBER_NOISEMODEL_H
#define __FABBER_NOISEMODEL_H 1

#include "dist_mvn.h"
#include "fwdmodel_linear.h"
#include "utils.h"

/**
 * Base class for parameters to a noise model
 *
 * Each derived NoiseModel will have a derived NoiseParams
 */
class NoiseParams
{
public:

	virtual NoiseParams* Clone() const = 0;
	virtual const NoiseParams& operator=(const NoiseParams& in) = 0;

	/**
	 * Output as a multivariate normal dist
	 */
	virtual const MVNDist OutputAsMVN() const = 0;

	virtual void InputFromMVN(const MVNDist& mvn) = 0;

	/**
	 * Dump human-readable debug output to default LOG
	 */
	virtual void Dump(const string indent = "") const = 0;

	virtual ~NoiseParams()
	{
	}
};

/**
 * Generic noise model
 *
 * This class & derived classes should be essentially data-free, instead storing
 * the relevant noise parameters in a NoiseParams-derived subclass
 *
 * FIXME unclear why params need to be separate?
 */
class NoiseModel
{

public:
	/**
	 * Static member function, to pick a noise model from a name
	 */
	static NoiseModel* NewFromName(const string& name);

	/**
	 * Create a new instance of this class. Subclasses
	 * implement this to produce an instance of themselves
	 *
	 * @return pointer to new instance.
	 */
	static NoiseModel* NewInstance();

	/**
	 * Initialize a new instance using configuration from the given
	 * arguments.
	 *
	 * @param args Configuration parameters.
	 */
	virtual void Initialize(FabberRunData& args)
	{
	}

	/**
	 * Create a new identical copy of this object (e.g. for spatial vb)
	 */
	//  virtual NoiseModel* Clone() const = 0;
	virtual NoiseParams* NewParams() const = 0;

	/**
	 * Load priors from file, and also initialize posteriors
	 */
	//  virtual void LoadPrior( const string& filename ) = 0;

	/**
	 * Suggest some nice default values for noise parameters:
	 */
	virtual void HardcodedInitialDists(NoiseParams& prior, NoiseParams& posterior) const = 0;

	/**
	 * Output your internal posterior distribution as an MVN.
	 * (bit of a hack -- some noise models don't fit into the MVN framework well)
	 */
	//  virtual const MVNDist GetResultsAsMVN() const = 0;

	/**
	 * Some noise models might want to precalculate things (for efficiency
	 * reasons), based on the length of the data... if you don't know what
	 * this is for then just ignore it.
	 */
	virtual void Precalculate(NoiseParams& noise, const NoiseParams& noisePrior, const NEWMAT::ColumnVector& sampleData) const
	{
	}

	virtual ~NoiseModel()
	{
	}

	// VB Updates

	// The following could potentially be split into substeps; but since
	// these would necessarily be model-specific, it's nice to have a
	// general catch-all update step.  Presumably this function
	// would call all the other functions in some order.

	virtual void UpdateNoise(NoiseParams& noise, const NoiseParams& noisePrior, const MVNDist& theta,
		const LinearFwdModel& model, const NEWMAT::ColumnVector& data) const = 0;

	virtual void UpdateTheta(const NoiseParams& noise,
			//    const NoiseParams& noisePrior,
			MVNDist& theta, const MVNDist& thetaPrior, const LinearFwdModel& model, const NEWMAT::ColumnVector& data,
			MVNDist* thetaWithoutPrior = NULL,
			// for --spatial-prior-output-correction
			float LMalpha = 0) const = 0;

	virtual double CalcFreeEnergy(const NoiseParams& noise, const NoiseParams& noisePrior, const MVNDist& theta,
		const MVNDist& thetaPrior, const LinearFwdModel& model, const NEWMAT::ColumnVector& data) const = 0;

	// ARD things
	double SetupARD(vector<int> ardindices, const MVNDist& theta, MVNDist& thetaPrior) const;

	double UpdateARD(vector<int> ardindices, const MVNDist& theta, MVNDist& thetaPrior) const;

	/**
	 * Return the number of noise parameters
	 */
	virtual int NumParams() = 0;

	// Potentially other functions could go here,
	// e.g. likelihood at a point (for MCMC) or sampling function (for Gibbs)

	//  virtual void SaveParams(const MVNDist& theta) { /* do nothing */ }
	//  virtual void RevertParams(MVNDist& theta)
	//        { throw Invalid_option("This noise model does not support reverting (don't use the trial-mode convergence detector with it)\n"); }

private:
	/**
	 * Prevent copying using anything other than the Clone() function.
	 * Could implement it, but not particularly useful and the default
	 * shallow copy is not right.
	 */
	const NoiseModel& operator=(const NoiseModel&) const
	{
		assert(false);
		return *this;
	}
};

/**
 * Handy mathematical function, used by some free energy calculations
 */
double gammaln(double xx);

/** 
 * \ref SingletonFactory that returns pointers to \ref NoiseModel.
 */
typedef SingletonFactory<NoiseModel> NoiseModelFactory;

#endif // __FABBER_NOISEMODEL_H
