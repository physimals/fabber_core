/*  fwdmodel_polynomial.h - Implements a simple polynomial model

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class PolynomialFwdModel: public FwdModel
{
public:
	static FwdModel* NewInstance();
	void Initialize(FabberRunData& args);
	virtual void Usage(std::ostream &stream) const;
	virtual void Evaluate(const ColumnVector& params, ColumnVector& result) const;

	virtual void DumpParameters(const ColumnVector& vec, const string& indents = "") const;
	virtual void NameParams(vector<string>& names) const;
	virtual int NumParams() const
	{
		return m_degree + 1;
	}
	string ModelVersion() const;
	void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

	virtual ~PolynomialFwdModel()
	{
	}
protected:
	int m_degree;
};

