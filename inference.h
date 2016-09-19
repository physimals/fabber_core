/*  inference.h - General inference technique base class

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#ifndef __FABBER_INFERENCE_H
#define __FABBER_INFERENCE_H 1

#include <map>
#include <string>
#include <typeinfo>
#include <vector>
#include <memory>

#include "dataset.h"
#include "fwdmodel.h"
#include "easylog.h"
#include "noisemodel.h"
#include "utils.h"

#ifdef __FABBER_MOTION
#include "Update_deformation.h"
#include "mcflirt/rigidreglib.h"
#endif //__FABBER_MOTION
class InferenceTechnique
{

public:
	/**
	 * Static member function, to pick an inference technique from a name
	 */
	static InferenceTechnique* NewFromName(const string& name);

	/**
	 * Create a new instance of this class.
	 * @return pointer to new instance.
	 */
	static InferenceTechnique* NewInstance();

	/**
	 * Constructor.
	 */
	InferenceTechnique() :
		model(NULL), noise(NULL)
	{
	}

	/**
	 * Initialize a new instance to use the given forward model
	 * and extract additional configuration from the given
	 * arguments.
	 * @param fwd_model Forward model to be used.
	 * @param args Additional configuration parameters.
	 */
	virtual void Initialize(FwdModel* fwd_model, FabberRunData& args);

	/**
	 * Perform inference using the given model upon the given data.
	 *
	 * This method should only be called after Initialize()
	 * Subclasses of InferenceTechnique must implement this method
	 * to carry out their given inference calculations
	 *
	 * @param data
	 */
	virtual void DoCalculations(FabberRunData& data) = 0;

	/**
	 * Save the results
	 */
	virtual void SaveResults(FabberRunData& data) const;

	/**
	 * Destructor.
	 */
	virtual ~InferenceTechnique();

protected:
	/**
	 * Pointer to forward model, passed in to initialize.
	 *
	 * Will not be deleted, that is the responsibility of
	 * the caller
	 */
	FwdModel* model;

	/**
	 * Number of model parameters.
	 *
	 * This is used regularly so it's sensible to keep a
	 * copy around
	 */
	int m_num_params;

	/**
	 * Number of noise parameters.
	 *
	 * This is used regularly so it's sensible to keep a
	 * copy around
	 */
	int m_noise_params;

	/**
	 * Noise model in use. This is created by the inference
	 * method deleted on destruction
	 */
	std::auto_ptr<NoiseModel> noise;

	bool saveModelFit;
	bool saveResiduals;

	/**
	 * Results of the inference method
	 *
	 * Vector of MVNDist, one for each voxel
	 * Each MVNDist contains the means and covariance/precisions for
	 * the parameters in the model
	 */
	vector<MVNDist*> resultMVNs;

	/**
	 * Used by Adrian's spatial priors research
	 */
	vector<MVNDist*> resultMVNsWithoutPrior;

	/** Free energy for each voxel? */
	vector<double> resultFs;

	void InitMVNFromFile(string continueFromFile, FabberRunData& allData, string paramFilename);

	/**
	 * Number of motion correction steps to run
	 */
	int Nmcstep;
private:
	/**
	 * Private to prevent assignment
	 */
	const InferenceTechnique& operator=(const InferenceTechnique& from)
	{
		assert(false);
		return from;
	}
};

/** 
 * \ref SingletonFactory that returns pointers to 
 * \ref InferenceTechnique.
 */
typedef SingletonFactory<InferenceTechnique> InferenceTechniqueFactory;

// Motion Correction class
#ifdef __FABBER_MOTION
//   NB: for now the mask should cover the *entire* image as we zero everything
//       outside of the mask, which is not good for registration
//       In future we'd need allData to be able to provide the original image (or something to)

class MCobj
{
public:
	MCobj(const FabberRunData& allData, int dof);
	void run_mc(const Matrix& modelpred_mat, Matrix& finalimage_mat);
	void set_num_iter(int nit)
	{	num_iter=nit;}
private:
	int userdof; // anything over 13 is full nonlinear
	int num_iter; // default 10
	volume<float> mask;
	Matrix affmat;
	mcflirt mcf;
	volume4D<float> defx;
	volume4D<float> defy;
	volume4D<float> defz;
	// things below are kept for efficiency (?) in order to avoid repeated allocation/destruction
	volume4D<float> tmpx;
	volume4D<float> tmpy;
	volume4D<float> tmpz;
	volume4D<float> modelpred;
	volume4D<float> finalimage;
	volume4D<float> wholeimage;
};

#endif // __FABBER_MOTION
#endif // __FABBER_INFERENCE_H
