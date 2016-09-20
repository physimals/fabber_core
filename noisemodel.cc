/*  noisemodel.cc - Class implementation for generic noise models

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include <math.h>

#include "noisemodel.h"

NoiseModel* NoiseModel::NewFromName(const string& name)
{
	NoiseModelFactory* factory = NoiseModelFactory::GetInstance();
	NoiseModel *noise = factory->Create(name);
	if (noise == NULL)
	{
		throw Invalid_option("Unrecognized --noise: " + name);
	}
	return noise;
}

// Calculate log-gamma from a Taylor expansion; good to one part in 2e-10.
double gammaln(double x)
{
	ColumnVector series(7);
	series << 2.5066282746310005 << 76.18009172947146 << -86.50532032941677 << 24.01409824083091 << -1.231739572450155
			<< 0.1208650973866179e-2 << -0.5395239384953e-5;

	double total = 1.000000000190015;
	for (int i = 2; i <= series.Nrows(); i++)
		total += series(i) / (x + i - 1);

	return log(series(1) * total / x) + (x + 0.5) * log(x + 5.5) - x - 5.5;
}

// ARD stuff
double NoiseModel::SetupARD(vector<int> ardindices, const MVNDist& theta, MVNDist& thetaPrior) const
{
	Tracer_Plus tr("Noisemodel::SetupARD");
	double Fard = 0;

	if (~ardindices.empty())
	{
		SymmetricMatrix PriorPrec;
		PriorPrec = thetaPrior.GetPrecisions();
		SymmetricMatrix PostCov = theta.GetCovariance();

		for (int i = 0; i < ardindices.size(); i++)
		{
			PriorPrec(ardindices[i], ardindices[i]) = 1e-12; //set prior to be initally non-informative
			thetaPrior.means(ardindices[i]) = 0;

			//set the Free energy contribution from ARD term
			double b = 2 / (theta.means(ardindices[i]) * theta.means(ardindices[i]) + PostCov(ardindices[i],
					ardindices[i]));
			Fard += -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5 * log(b); //taking c as 0.5 - which it will be!
		}

		thetaPrior.SetPrecisions(PriorPrec);
	}

	return Fard;
}

double NoiseModel::UpdateARD(vector<int> ardindices, const MVNDist& theta, MVNDist& thetaPrior) const
{
	Tracer_Plus tr("Noisemodel::UpdateARD");
	double Fard = 0;

	if (~ardindices.empty())
	{
		SymmetricMatrix PriorCov;
		SymmetricMatrix PostCov;
		PriorCov = thetaPrior.GetCovariance();
		PostCov = theta.GetCovariance();

		for (int i = 0; i < ardindices.size(); i++)
		{
			PriorCov(ardindices[i], ardindices[i]) = theta.means(ardindices[i]) * theta.means(ardindices[i]) + PostCov(
					ardindices[i], ardindices[i]);

			//set the Free energy contribution from ARD term
			double b = 2 / (theta.means(ardindices[i]) * theta.means(ardindices[i]) + PostCov(ardindices[i],
					ardindices[i]));
			Fard += -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5 * log(b); //taking c as 0.5 - which it will be!
		}

		thetaPrior.SetCovariance(PriorCov);
	}

	return Fard;
}
