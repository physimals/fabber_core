/*  noisemodel_ar.cc - Class implementation for the AR(1) noise model

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */
#pragma once

#include "noisemodel.h"

#include "dist_gamma.h"

#include <ostream>
#include <string>
#include <vector>

class Ar1cParams;

// Helper class -- caches some of the AR matrices
class Ar1cMatrixCache
{
public:
    explicit Ar1cMatrixCache(int numPhis);
    Ar1cMatrixCache(const Ar1cMatrixCache &from);
    const NEWMAT::SymmetricBandMatrix &GetMatrix(unsigned n, unsigned a12pow, unsigned a3pow) const;
    const NEWMAT::SymmetricBandMatrix &GetMarginal(unsigned n) const;
    void Update(const Ar1cParams &dist, int nTimes);

private:
    std::vector<NEWMAT::SymmetricBandMatrix> alphaMarginals;
    // recalculated whenever alpha changes
    unsigned FlattenIndex(unsigned n, unsigned a12pow, unsigned a34pow) const;

    // should only be calculated once
    // Note that if more than one model is being inferred upon at a time,
    // this will be unnecessarily duplicated in every one of them --
    // might speed things up considerably by sharing.
    std::vector<NEWMAT::SymmetricBandMatrix> alphaMatrices;

    int nPhis;
};

// Parameter-storage class -- it's really just an enhanced structure
class Ar1cParams : public NoiseParams
{
public:
    // Constructors
    Ar1cParams(int nAlpha, int nPhi);
    Ar1cParams(const Ar1cParams &from);
    virtual const Ar1cParams &operator=(const NoiseParams &in);
    virtual Ar1cParams *Clone() const;

    virtual const MVNDist OutputAsMVN() const;
    virtual void InputFromMVN(const MVNDist &mvn);
    virtual void Dump(std::ostream &os) const;

private:
    friend class Ar1cNoiseModel; // Needs to use this class like it's a structure
    friend class Ar1cMatrixCache;
    MVNDist alpha;
    std::vector<GammaDist> phis;

    Ar1cMatrixCache alphaMat;
};

class Ar1cNoiseModel : public NoiseModel
{
public:
    static NoiseModel *NewInstance();

    virtual void Initialize(FabberRunData &args);
    virtual Ar1cParams *NewParams() const;
    int NumParams();

    virtual void HardcodedInitialDists(NoiseParams &prior, NoiseParams &posterior) const;

    /** Used to pre-evaluate the alpha matrices in the cache */
    virtual void Precalculate(NoiseParams &noise, const NoiseParams &noisePrior,
        const NEWMAT::ColumnVector &sampleData) const;

    virtual void UpdateNoise(NoiseParams &noise, const NoiseParams &noisePrior,
        const MVNDist &theta, const LinearFwdModel &linear, const NEWMAT::ColumnVector &data) const;

    virtual void UpdateTheta(const NoiseParams &noise, MVNDist &theta, const MVNDist &thetaPrior,
        const LinearFwdModel &model, const NEWMAT::ColumnVector &data,
        MVNDist *thetaWithoutPrior = NULL, float LMalpha = 0) const;

    virtual double CalcFreeEnergy(const NoiseParams &noise, const NoiseParams &noisePrior,
        const MVNDist &theta, const MVNDist &thetaPrior, const LinearFwdModel &model,
        const NEWMAT::ColumnVector &data) const;

protected:
    std::string ar1Type;
    int NumAlphas() const; // converts the above string into a number
    int nPhis;

    virtual void UpdateAlpha(NoiseParams &noise, const NoiseParams &noisePrior,
        const MVNDist &theta, const LinearFwdModel &model, const NEWMAT::ColumnVector &data) const;

    virtual void UpdatePhi(NoiseParams &noise, const NoiseParams &noisePrior, const MVNDist &theta,
        const LinearFwdModel &model, const NEWMAT::ColumnVector &data) const;
};
