/*  fwdmodel_polynomial.h - Implements a simple polynomial model

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>

/**
 * Forward model which fits to a simple polynomial function
 * e.g. a + bx + cx^2 + ....
 *
 * The parameter 'degree' can be used to set the maximum power
 * in the polynomial, e.g. if degree=3, there will be 4 parameters
 */
class PolynomialFwdModel: public FwdModel
{
public:
	static FwdModel* NewInstance();
	void Initialize(FabberRunData& args);
	virtual std::vector<OptionSpec> GetOptions() const;

	virtual void Evaluate(const ColumnVector& params, ColumnVector& result) const;

	virtual void DumpParameters(const ColumnVector& vec, const std::string& indents = "") const;
	virtual void NameParams(std::vector<std::string>& names) const;
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

