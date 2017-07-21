#include "transforms.h"

DistParams Transform::ToModel(DistParams params) const
{
    return DistParams(ToModel(params.mean()), ToModel(params.var()));
}

DistParams Transform::ToFabber(DistParams params) const
{
    return DistParams(ToFabber(params.mean()), ToFabber(params.var()));
}

DistParams LogTransform::ToModel(DistParams params) const
{
    double mean = exp(params.mean() + params.var()/2);
    double var = (exp(params.var())-1)*exp(2*params.mean() + params.var());
    return DistParams(mean, var);
}

DistParams LogTransform::ToFabber(DistParams params) const
{
    double mean = 2*log(params.mean()) - 0.5*log(params.var() + params.mean()*params.mean());
    double var = log(params.var()/(params.mean()*params.mean()) + 1);
    return DistParams(mean, var);
}

DistParams FractionalTransform::ToModel(DistParams params) const
{
    double mean = exp(params.mean());
    double var = exp(params.var());
    return DistParams(mean, var);
}

DistParams FractionalTransform::ToFabber(DistParams params) const
{
    double mean = log(params.mean());
    double var = log(params.var());
    return DistParams(mean, var);
}

/** Singleton instance of identity transform */
const Transform *TRANSFORM_IDENTITY()
{
    static IdentityTransform t;
    return &t;
}

/** Singleton instance of log transform */
const Transform *TRANSFORM_LOG()
{
    static LogTransform t;
    return &t;
}

/** Singleton instance of softplus transform */
const Transform *TRANSFORM_SOFTPLUS()
{
    static SoftPlusTransform t;
    return &t;
}

/** Singleton instance of fractional transform */
const Transform *TRANSFORM_FRACTIONAL()
{
    static FractionalTransform t;
    return &t;
}

const Transform *GetTransform(std::string id)
{
    if (id == TRANSFORM_CODE_IDENTITY) return TRANSFORM_IDENTITY();
    else if (id == TRANSFORM_CODE_LOG) return TRANSFORM_LOG();
    else if (id == TRANSFORM_CODE_SOFTPLUS) return TRANSFORM_SOFTPLUS();
    else if (id == TRANSFORM_CODE_FRACTIONAL) return TRANSFORM_FRACTIONAL();
    else throw InvalidOptionValue("PSP_byname<n>_transform", id, "Supported transforms: I, L, S");
}
