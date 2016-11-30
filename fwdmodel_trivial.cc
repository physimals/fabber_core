/*  fwdmodel_trivial.cc - Implements an utterly uninteresting model

 Martin Craig

 Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_trivial.h"

using NEWMAT::ColumnVector;
using namespace std;

FwdModel* TrivialFwdModel::NewInstance()
{
	return new TrivialFwdModel();
}

std::string TrivialFwdModel::GetDescription() const
{
	return "Trivial forward model which fits to a constant function. Mainly useful for testing";
}

string TrivialFwdModel::ModelVersion() const
{
	return "1.0";
}

void TrivialFwdModel::Initialize(FabberRunData& args)
{
}

int TrivialFwdModel::NumParams() const
{
	return 1;
}

void TrivialFwdModel::NameParams(vector<string>& names) const
{
	names.clear();
	names.push_back("p");
}


void TrivialFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
	// output = parameter
	assert(params.Nrows() == 1);
	result.ReSize(data.Nrows());
	result = params(1);
}

void TrivialFwdModel::HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const
{
	// Have to implement this
	assert(prior.means.Nrows() == 1);
	prior.means = 0;
	prior.SetPrecisions(NEWMAT::IdentityMatrix(1) * 1e-12);
	posterior = prior;
}
