/*  noisemodel.cc - Class implementation for generic noise models

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "noisemodel.h"

#include "rundata.h"
#include "tools.h"

#include <miscmaths/miscmaths.h>

#include <math.h>
#include <string>

using MISCMATHS::digamma;
using NEWMAT::SymmetricMatrix;
using namespace std;

NoiseModel *NoiseModel::NewFromName(const string &name)
{
    NoiseModelFactory *factory = NoiseModelFactory::GetInstance();
    NoiseModel *noise = factory->Create(name);
    if (noise == NULL)
    {
        throw InvalidOptionValue("noise", name, "Unrecognized noise type");
    }
    return noise;
}

void NoiseModel::Initialize(FabberRunData &rundata)
{
    m_log = rundata.GetLogger();

    // Read masked time points option if any have been specified
    m_masked_tpoints = rundata.GetIntList("mt", 1);
}
// ARD stuff
double NoiseModel::SetupARD(vector<int> ardindices, const MVNDist &theta, MVNDist &thetaPrior) const
{
    double Fard = 0;

    if (!ardindices.empty())
    {
        SymmetricMatrix PriorPrec;
        PriorPrec = thetaPrior.GetPrecisions();
        SymmetricMatrix PostCov = theta.GetCovariance();

        for (size_t i = 0; i < ardindices.size(); i++)
        {
            PriorPrec(ardindices[i], ardindices[i])
                = 1e-12; // set prior to be initally non-informative
            thetaPrior.means(ardindices[i]) = 0;

            // set the Free energy contribution from ARD term
            double b = 2 / (theta.means(ardindices[i]) * theta.means(ardindices[i])
                               + PostCov(ardindices[i], ardindices[i]));
            Fard += -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5)
                - 0.5 * log(b); // taking c as 0.5 - which it will be!
        }

        thetaPrior.SetPrecisions(PriorPrec);
    }

    return Fard;
}

double NoiseModel::UpdateARD(
    vector<int> ardindices, const MVNDist &theta, MVNDist &thetaPrior) const
{
    double Fard = 0;

    if (!ardindices.empty())
    {
        SymmetricMatrix PriorCov;
        SymmetricMatrix PostCov;
        PriorCov = thetaPrior.GetCovariance();
        PostCov = theta.GetCovariance();

        for (size_t i = 0; i < ardindices.size(); i++)
        {
            PriorCov(ardindices[i], ardindices[i])
                = theta.means(ardindices[i]) * theta.means(ardindices[i])
                + PostCov(ardindices[i], ardindices[i]);

            // set the Free energy contribution from ARD term
            double b = 2 / (theta.means(ardindices[i]) * theta.means(ardindices[i])
                               + PostCov(ardindices[i], ardindices[i]));
            Fard += -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5)
                - 0.5 * log(b); // taking c as 0.5 - which it will be!
        }

        thetaPrior.SetCovariance(PriorCov);
    }

    return Fard;
}
