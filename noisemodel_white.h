/*  noisemodel_ar.h - Class declaration for the multiple white noise model

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "noisemodel.h"

#include "dist_gamma.h"

/**
 * Parameters for the white noise model
 *
 * The noise model can have any number of parameters (Phis)
 * which apply to different samples in the timeseries, e.g.
 * a noise pattern of 12121212... has two parameters
 * one applying to odd numbered samples, and one to even.
 *
 * Each Phi is associated with a gamma distribution
 */
class WhiteParams: public NoiseParams
{
public:
	virtual WhiteParams* Clone() const
	{
		return new WhiteParams(*this);
	}

	virtual const WhiteParams& operator=(const NoiseParams& in)
	{
		const WhiteParams& from = dynamic_cast<const WhiteParams&> (in);
		assert(nPhis == from.nPhis);
		phis = from.phis;
		return *this;
	}

	virtual const MVNDist OutputAsMVN() const;
	virtual void InputFromMVN(const MVNDist& mvn);

	virtual void Dump(const string indent = "") const;

	WhiteParams(int N) :
		nPhis(N), phis(N)
	{
		return;
	}
	WhiteParams(const WhiteParams& from) :
		nPhis(from.nPhis), phis(from.phis)
	{
		return;
	}

private:
	friend class WhiteNoiseModel;
	const int nPhis;
	vector<GammaDist> phis;
};

class WhiteNoiseModel: public NoiseModel
{
public:
	/**
	 * Create a new instance of white noise. Used by factory
	 * to create noise models by name. Parameters are set
	 * during initialization
	 */
	static NoiseModel* NewInstance();

	virtual WhiteParams* NewParams() const
	{
		return new WhiteParams(Qis.size());
	}

	virtual void HardcodedInitialDists(NoiseParams& prior, NoiseParams& posterior) const;

	virtual void Initialize(FabberRunData& args);

	virtual ~WhiteNoiseModel()
	{
	}

	/**
	 * Update the noise parameters
	 */
	virtual void UpdateNoise(NoiseParams& noise, const NoiseParams& noisePrior, const MVNDist& theta,
		const LinearFwdModel& model, const NEWMAT::ColumnVector& data) const;

	/**
	 * Update the model parameters?
	 */
	virtual void
	UpdateTheta(const NoiseParams& noise, MVNDist& theta, const MVNDist& thetaPrior, const LinearFwdModel& model,
	const NEWMAT::ColumnVector& data, MVNDist* thetaWithoutPrior = NULL, float LMalpha = 0) const;

	virtual double CalcFreeEnergy(const NoiseParams& noise, const NoiseParams& noisePrior, const MVNDist& theta,
		const MVNDist& thetaPrior, const LinearFwdModel& model, const NEWMAT::ColumnVector& data) const;

	int NumParams();

protected:
	string phiPattern;

	double lockedNoiseStdev; // Allow phi to be locked externally
	double phiprior; //allow external setting of the prior nosie std deviation (and thence phi)

	// Diagonal matrices, indicating which data points use each phi
	mutable vector<NEWMAT::DiagonalMatrix> Qis; // mutable because it's used as a cache
	void MakeQis(int dataLen) const;
};
