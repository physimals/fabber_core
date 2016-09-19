/*  fwdmodel_trivial.cc - Implements a polynomial model

 Martin Craig

 Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include <iostream>

#include "fwdmodel_poly.h"

void PolynomialFwdModel::Usage(std::ostream &stream) const
{
	stream << "--degree=<polynomial degree> - Required" << endl;
}

FwdModel* PolynomialFwdModel::NewInstance()
{
	return new PolynomialFwdModel();
}

void PolynomialFwdModel::Initialize(FabberRunData& args)
{
	m_degree = convertTo<int> (args.GetString("degree"));
	//	m_coeffs.resize(m_degree);
}

string PolynomialFwdModel::ModelVersion() const
{
	return "0"; // This model does not really deserve a version
}

void PolynomialFwdModel::Evaluate(const ColumnVector& params,
		ColumnVector& result) const
{
	assert(params.Nrows() == m_degree+1);
	result.ReSize(data.Nrows());

	for (int i = 1; i <= data.Nrows(); i++)
	{
		double res = 0;
		int p = 1;
		for (int n = 0; n <= m_degree; n++)
		{
			res += params(n+1) * p;
			p *= i;
		}
		result(i) = res;
	}
}

void PolynomialFwdModel::DumpParameters(const ColumnVector& vec,
		const string& indent) const
{
	LOG << indent << vec << endl;
}

void PolynomialFwdModel::HardcodedInitialDists(MVNDist& prior,
		MVNDist& posterior) const
{
	// Have to implement this
	assert(prior.means.Nrows() == m_degree+1);
	prior.means = 0;
	prior.SetPrecisions(IdentityMatrix(m_degree + 1) * 1e-12);
	posterior = prior;
}

void PolynomialFwdModel::NameParams(vector<string>& names) const
{
	names.clear();
	for (int i = 0; i <= m_degree; i++)
	{
		names.push_back((string)"c"+stringify(i));
	}
}

