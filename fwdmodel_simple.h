/*  fwdmodel_simple.h - Implements the simplified ASL model

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"
#include "inference.h"
#include <string>
using namespace std;

class SimpleFwdModelIdStruct
{
public:
	SimpleFwdModelIdStruct()
	{
		Nbasis = 0;
	}
	void Define(int NumBasisFcns)
	{
		assert(Nbasis==0);
		assert(NumBasisFcns>0);
		Nbasis = NumBasisFcns;
	}

	//  virtual void DumpVector(const ColumnVector& vec, const string& indent = "");
	//  virtual ~SimpleFwdModelIdStruct() { return; }

	// N=3: 1=Q0 2=Q1 3=Q2 4=Q3 5=U0 6=U1 7=U2 8=U3 9=R0 10=R1 11=R2 12=R3
	double Q0(const ColumnVector& p) const
	{
		return p(1);
	}
	ReturnMatrix Qn(const ColumnVector& p) const
	{
		AssertValid(p);
		ColumnVector tmp = p.Rows(2, Nbasis + 1);
		return tmp;
	}
	double U0(const ColumnVector& p) const
	{
		return p(2 + Nbasis);
	}
	ReturnMatrix Un(const ColumnVector& p) const
	{
		AssertValid(p);
		ColumnVector tmp = p.Rows(3 + Nbasis, 2 + 2 * Nbasis);
		return tmp;
	}
	double R0(const ColumnVector& p) const
	{
		AssertValid(p);
		return p(3 + 2 * Nbasis);
	}
	ReturnMatrix Rn(const ColumnVector& p) const
	{
		AssertValid(p);
		ColumnVector tmp = p.Rows(4 + 2 * Nbasis, 3 + 3 * Nbasis);
		return tmp;
	}
	void NameParams(vector<string>& names) const;

private:
	int Nbasis;
	void AssertValid(const ColumnVector& p) const
	{
		assert(Nbasis>0);
		assert(p.Nrows() == 3 + 3*Nbasis);
	}
};

class SimpleFwdModel: public FwdModel
{
public:
	// Virtual function overrides
	virtual void Evaluate(const ColumnVector& params, ColumnVector& result) const;

	virtual void DumpParameters(const ColumnVector& vec, const string& indents = "") const;
	virtual void NameParams(vector<string>& names) const
	{
		id.NameParams(names);
	}
	virtual int NumParams() const
	{
		return 3 + 3 * basis.Ncols();
	}
	string ModelVersion() const;

	virtual ~SimpleFwdModel()
	{
		return;
	}

	virtual void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const
	{
		prior.means = 0;
		prior.SetPrecisions(IdentityMatrix(NumParams()) * 1e-12);
		posterior = prior;
		// Informative starting points
		// I'm not very happy with this whole setup... the physical arrangment
		// of elements in the vector should only appear once.
		posterior.means(1) = 50;
		posterior.means(2 + basis.Ncols()) = 1.5e4;
		posterior.means(3 + 2 * basis.Ncols()) = 25;
		assert(id.Q0(posterior.means) == 50);
		assert(id.U0(posterior.means) == 1.5e4);
		assert(id.R0(posterior.means) == 25);
	}

	// Constructor
	SimpleFwdModel(ArgsType& args);

	//  SimpleFwdModel(const SimpleFwdModel& from); // copy constructor - default ok?

protected:
	// Constants
	ColumnVector echoTime;
	Matrix basis; // basis in columns?
	ColumnVector rho;
	SimpleFwdModelIdStruct id;
};

