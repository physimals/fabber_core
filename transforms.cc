#include "transforms.h"

DistParams Transform::ToModel(DistParams params) const
{
    double mean = ToModel(params.mean());
    double var = ToModelVar(params.var());
    return DistParams(mean, var);
}

DistParams Transform::ToFabber(DistParams params) const
{
    double mean = ToFabber(params.mean());
    double var = ToFabberVar(params.var());
    return DistParams(mean, var);
}

double Transform::ToModelVar(double val) const
{
    return pow(ToModel(sqrt(val)) - ToModel(0), 2);
}

double Transform::ToFabberVar(double val) const
{
    return pow(ToFabber(ToModel(0) + sqrt(val)), 2);
}

const Transform *TRANSFORM_IDENTITY()
{
    static IdentityTransform t;
    return &t;
}

const Transform *TRANSFORM_LOG()
{
    static LogTransform t;
    return &t;
}

const Transform *TRANSFORM_SOFTPLUS()
{
    static SoftPlusTransform t;
    return &t;
}

const Transform *TRANSFORM_FRACTIONAL()
{
    static FractionalTransform t;
    return &t;
}

const Transform *TRANSFORM_ABS()
{
    static AbsTransform t;
    return &t;
}

const Transform *GetTransform(std::string id)
{
    if (id == TRANSFORM_CODE_IDENTITY)
        return TRANSFORM_IDENTITY();
    else if (id == TRANSFORM_CODE_LOG)
        return TRANSFORM_LOG();
    else if (id == TRANSFORM_CODE_SOFTPLUS)
        return TRANSFORM_SOFTPLUS();
    else if (id == TRANSFORM_CODE_FRACTIONAL)
        return TRANSFORM_FRACTIONAL();
    else if (id == TRANSFORM_CODE_ABS)
        return TRANSFORM_ABS();
    else
        throw InvalidOptionValue("PSP_byname<n>_transform", id, "Supported transforms: I, L, S, F");
}
