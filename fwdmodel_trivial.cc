/*  fwdmodel_trivial.cc - Implements an utterly uninteresting model

 Martin Craig

 Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include <iostream>

#include "fwdmodel_trivial.h"

void TrivialFwdModel::Usage(std::ostream &stream) const
{
    stream << "No configurable parameters at all!" << endl;
}

FwdModel* TrivialFwdModel::NewInstance()
{
    return new TrivialFwdModel();
}

void TrivialFwdModel::Initialize(FabberRunData& args)
{
}

string TrivialFwdModel::ModelVersion() const
{
    return "0"; // This model does not really deserve a version
}

void TrivialFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
    // output = parameter
    assert(params.Nrows() == 1);
    result.ReSize(data.Nrows());
    result = params(1);
}

void TrivialFwdModel::DumpParameters(const ColumnVector& vec, const string& indent) const
{
     LOG << indent << vec << endl;
}

void TrivialFwdModel::HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const
{
    // Have to implement this
    assert(prior.means.Nrows() == 1);
    prior.means = 0;
    prior.SetPrecisions(IdentityMatrix(1) * 1e-12);
    posterior = prior;
}

void TrivialFwdModel::NameParams(vector<string>& names) const
{
    names.clear();
    names.push_back("p");
}

