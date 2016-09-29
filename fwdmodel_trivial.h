/*  fwdmodel_simple.h - Implements the simplified ASL model

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"

#include <string>
#include <vector>

class TrivialFwdModel: public FwdModel
{
public:
	static FwdModel* NewInstance();
	std::string GetDescription() const;
	string ModelVersion() const;

	void Initialize(FabberRunData& args);
	void NameParams(std::vector<std::string>& names) const;
	int NumParams() const;

	void Evaluate(const NEWMAT::ColumnVector& params, NEWMAT::ColumnVector& result) const;
	void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;
};

