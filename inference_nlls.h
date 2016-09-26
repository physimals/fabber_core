/* inference_nlls.h - Non-Linear Least Squares class declarations

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford */

/*  CCOPYRIGHT  */

#include "inference.h"
#include "miscmaths/nonlin.h"
#include <boost/shared_ptr.hpp>
#include "miscmaths/bfmatrix.h"

/**
 * Inference technique using non-linear least squares
 */
class NLLSInferenceTechnique: public InferenceTechnique
{
public:
	/**
	 * Create a new NLLSInferenceTechnique instance
	 *
	 * This infers model parameters by minimising the sum of
	 * the squares of the difference between the data and
	 * the model for each sample in the timeseries.
	 *
	 * This calculation is done independently for each voxel
	 */
	static InferenceTechnique* NewInstance();

	virtual vector<OptionSpec> GetOptions() const;
	virtual std::string GetDescription() const;
	virtual string GetVersion() const;

	virtual void Initialize(FwdModel* fwd_model, FabberRunData& args);
	virtual void DoCalculations(FabberRunData& data);
	virtual ~NLLSInferenceTechnique();
protected:
	const MVNDist* initialFwdPosterior;
	bool vbinit;
	bool lm;
};

/**
 * Cost function for NLLS method
 */
class NLLSCF: public NonlinCF
{
public:
	/**
	 * Create a cost function evaluator
	 *
	 * @param pdata Data we are trying to fit to
	 * @param pm  Model whose parameters we are trying to
	 *            adjust to fit the data
	 */
	NLLSCF(const ColumnVector& pdata, const FwdModel* pm) :
		y(pdata), model(pm), linear(pm)
	{
	}

	/**
	 * Virtual destructor in case subclasses are created
	 */
	virtual ~NLLSCF()
	{
	}

	/**
	 * Calculate the cost function
	 *
	 * This is the sum of the squares of the differences
	 * between the model and the data. The model is
	 * evalulated for the current parameters and compared
	 * to the data provided in the cost function constructor.
	 *
	 * @param p Current parameters
	 */
	virtual double cf(const ColumnVector& p) const;

	/**
	 * Calculate the gradient of the linearized model
	 *
	 * @param p Current model parameters
	 */
	virtual ReturnMatrix grad(const ColumnVector& p) const;

	/**
	 * Calculate the Hessian of the linearized model
	 *
	 * @param p Current model parameters
	 * @param iptr FIXME ???
	 */
	virtual boost::shared_ptr<BFMatrix> hess(const ColumnVector& p, boost::shared_ptr<BFMatrix> iptr) const;
private:
	const ColumnVector y;
	const FwdModel* model;
	mutable LinearizedFwdModel linear;
};
