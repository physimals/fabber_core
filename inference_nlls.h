#ifndef NO_NLLS
/* inference_nlls.h - Non-Linear Least Squares class declarations

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford */

/*  CCOPYRIGHT  */
#pragma once

#include "inference.h"

#include <miscmaths/bfmatrix.h>
#include <miscmaths/nonlin.h>

#include <boost/shared_ptr.hpp>

/**
 * Inference technique using non-linear least squares
 */
class NLLSInferenceTechnique : public InferenceTechnique
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
    static InferenceTechnique *NewInstance();

    virtual void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;
    virtual std::string GetVersion() const;

    virtual void Initialize(FwdModel *fwd_model, FabberRunData &args);
    virtual void DoCalculations(FabberRunData &data);

protected:
    const MVNDist *initialFwdPosterior;
    bool m_vbinit;
    bool m_lm;
};

/**
 * Cost function for NLLS method
 */
class NLLSCF : public MISCMATHS::NonlinCF
{
public:
    /**
     * Create a cost function evaluator
     *
     * @param pdata Data we are trying to fit to
     * @param pm  Model whose parameters we are trying to
     *            adjust to fit the data
     */
    NLLSCF(const NEWMAT::ColumnVector &pdata, const FwdModel *pm)
        : m_data(pdata)
        , m_model(pm)
        , m_linear(pm)
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
    virtual double cf(const NEWMAT::ColumnVector &p) const;

    /**
     * Calculate the gradient of the linearized model
     *
     * @param p Current model parameters
     */
    virtual NEWMAT::ReturnMatrix grad(const NEWMAT::ColumnVector &p) const;

    /**
     * Calculate the Hessian of the linearized model
     *
     * @param p Current model parameters
     * @param iptr FIXME ???
     */
    virtual boost::shared_ptr<MISCMATHS::BFMatrix> hess(
        const NEWMAT::ColumnVector &p, boost::shared_ptr<MISCMATHS::BFMatrix> iptr) const;

private:
    const NEWMAT::ColumnVector m_data;
    const FwdModel *m_model;
    mutable LinearizedFwdModel m_linear;
};
#endif
