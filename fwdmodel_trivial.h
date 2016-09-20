/*  fwdmodel_simple.h - Implements the simplified ASL model

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class TrivialFwdModel: public FwdModel
{
public:
	static FwdModel* NewInstance();
	void Initialize(FabberRunData& args);
	virtual void Evaluate(const ColumnVector& params, ColumnVector& result) const;

	virtual void DumpParameters(const ColumnVector& vec, const string& indents = "") const;
	virtual void NameParams(vector<string>& names) const;
	virtual int NumParams() const
	{
		return 1;
	}
	string ModelVersion() const;
	void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;

	virtual ~TrivialFwdModel()
	{
		return;
	}

};

