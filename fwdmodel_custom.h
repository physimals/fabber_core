/*  fwdmodel_custom.h - A place for your own quick'n'dirty forward model implementations.

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT */

/* fwdmodel_custom.h
 * Class declaration for a custom forward model
 * Template written by Adrian Groves, 2008
 * FMRIB Centre, University of Oxford
 *
 * Last modified: $Date: 2009/01/23 12:58:53 $ $Author: adriang $ $Revision: 1.1 $
 */

#ifndef __FABBER_FWDMODEL_CUSTOM_H
#define __FABBER_FWDMODEL_CUSTOM_H 1

#include "fwdmodel.h"

using namespace NEWMAT;
using namespace std;

/* Ideally you should create your own fwdmodel_*.h and .cc files to implement your models.
 * This "custom" model gives you a place to put your code and quickly start fitting your model.
 */

class CustomFwdModel: public FwdModel
{
public:
	CustomFwdModel(ArgsType& args);
	virtual ~CustomFwdModel()
	{
		return;
	}

	virtual void Evaluate(const ColumnVector& params, ColumnVector& result) const;
	virtual string ModelVersion() const;
	virtual int NumParams() const;
	virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;
	virtual void NameParams(vector<string>& names) const;
	// use default implementation of DumpParameters

	static void Usage(std::ostream &stream);

protected:
	// TODO: put your implementation-specific variables here.  These include any constants or
	// flags that are read from the command line.

	// int numRepeats;          // e.g. number of repeats
	// double T1;               // e.g. a fixed scan parameter
	// ColumnVector echoTimes;  // e.g. a list of TEs
	// Matrix basis;            // e.g. a basis set
	// bool wantFriesWithThat;  // e.g. a boolean option
};

#endif // __FABBER_FWDMODEL_CUSTOM_H
