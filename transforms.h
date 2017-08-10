#pragma once
/**
 * transforms.h
 *
 * Classes defining transformations on model parameters. These may be used
 * as a means of preventing model parameters from violating hard limits,
 * such as not becoming negative. Transformations are assumed to be
 * stateless, hence singleton instances can be used.
 *
 * Copyright (C) 2007-2017 University of Oxford
 */

/*  CCOPYRIGHT */

#include "rundata.h"

#include <math.h>
#include <string>

const std::string TRANSFORM_CODE_IDENTITY = "I";
const std::string TRANSFORM_CODE_LOG = "L";
const std::string TRANSFORM_CODE_SOFTPLUS = "S";
const std::string TRANSFORM_CODE_FRACTIONAL = "F";

/**
 * Immutable object describing distribution parameters for single model parameter
 */
struct DistParams
{
    DistParams(double m = 0, double v = 1)
        : m_mean(m)
        , m_var(v)
        , m_prec(1 / v)
    {
    }
    double mean() const { return m_mean; }
    double var() const { return m_var; }
    double prec() const { return m_prec; }
private:
    double m_mean;
    double m_var;
    double m_prec;
};

/**
 * Abstract base class for parameter transformations.
 */
class Transform
{
public:
    virtual ~Transform() {}
    /**
     * Transform the Fabber internal value (which is assumed to have a Gaussian
     * distribution) to the value required by the model
     */
    virtual double ToModel(double val) const = 0;

    /**
     * Transform a model parameter value to the value to be used by Fabber internally
     * (and modelled as a Gaussian random variable)
     */
    virtual double ToFabber(double val) const = 0;

    /**
     * Transform the Fabber internal mean/variance (of the Gaussian
     * distribution) to the mean/variance of the model distribution
     *
     * The default just applies ToModel to the mean and performs
     * an ad-hoc transformation to the variance by taking the square
     * root, adding it to the mean, transforming the result, then
     * subtracting the transformed mean and squaring. This just ensures that
     * mean + s.d. transforms into transformed mean + transformed s.d.
     *
     * Particular transformations may implement this differently, e.g.
     * the log transformation can use the known expression for the
     * variance of the log-normal distribution in terms of the
     * underlying normal distribution.
     */
    virtual DistParams ToModel(DistParams params) const;

    /**
     * Transform the mean of the model distribution into the
     * mean of the corresponding Gaussian distribution.
     *
     * The default just applies ToFabber to the mean and transforms
     * the variance in an ad-hoc manner - @see ToModel for details.
     */
    virtual DistParams ToFabber(DistParams params) const;
};

/**
 * Trivial implementation which performs no transformation
 */
class IdentityTransform : public Transform
{
public:
    double ToModel(double val) const { return val; }
    double ToFabber(double val) const { return val; }
    DistParams ToModel(DistParams params) const { return params; }
    DistParams ToFabber(DistParams params) const { return params; }
};

/**
 * Transformation for a parameter which is log-normal
 */
class LogTransform : public Transform
{
public:
    virtual DistParams ToModel(DistParams params) const;
    virtual DistParams ToFabber(DistParams params) const;
    double ToModel(double val) const { return exp(val); }
    double ToFabber(double val) const { return log(val); }
};

/**
 * 'Softplus' transformation
 *
 * This is an alternative to the log transform for a parameter which
 * is strictly positive. For positive values it approaches the
 * identity transformation so avoiding any issues associated with the
 * rapid growth of the exponential function.
 */
class SoftPlusTransform : public Transform
{
public:
    double ToModel(double val) const { return log(1 + exp(val)); }
    double ToFabber(double val) const { return log(exp(val) - 1); }
};

/**
 * 'Fractional' transformation
 *
 * This enforces parameter values between 0 and 1, exclusive
 * The variance is not transformed at all. This transformation
 * is intended for parameters which are not biophysical and
 * therefore we are not interested in specifying informative
 * priors, or analysing output variance.
 */
class FractionalTransform : public Transform
{
public:
    DistParams ToModel(DistParams params) const;
    DistParams ToFabber(DistParams params) const;
    double ToModel(double val) const { return 1 / (1 + exp(val)); }
    double ToFabber(double val) const { return log(1 / val - 1); }
};

/** Singleton instance of identity transform */
const Transform *TRANSFORM_IDENTITY();

/** Singleton instance of log transform */
const Transform *TRANSFORM_LOG();

/** Singleton instance of softplus transform */
const Transform *TRANSFORM_SOFTPLUS();

/** Singleton instance of fractional transform */
const Transform *TRANSFORM_FRACTIONAL();

const Transform *GetTransform(std::string id);
