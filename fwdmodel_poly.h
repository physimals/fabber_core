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
 *
 * Note that this class is mostly for testing purposes and is not
 * designed to be overriden.
 */
class PolynomialFwdModel: public FwdModel
{
public:
	static FwdModel* NewInstance();
	void GetOptions(std::vector<OptionSpec> &opts) const;
	std::string GetDescription() const;
	std::string ModelVersion() const;

	void Initialize(FabberRunData& args);
	int NumParams() const;
	void NameParams(std::vector<std::string>& names) const;

	void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;
	void Evaluate(const ColumnVector& params, ColumnVector& result) const;

private:
	int m_degree;
};

