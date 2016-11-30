/*  fwdmodel_custom.cc - A place for your own quick'n'dirty forward model implementations.

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT */

/* fwdmodel_custom.h
 * Implementation for a custom forward model class
 * Template written by Adrian Groves, 2008
 * FMRIB Centre, University of Oxford
 *
 * Last modified: $Date: 2012/01/13 12:00:59 $ $Author: adriang $ $Revision: 1.2 $
 */

#include "fwdmodel_custom.h"
#include "utils/tracer_plus.h"

using Utilities::Tracer_Plus;

// Constructor:
CustomFwdModel::CustomFwdModel(ArgsType& args)
{
	// Read any command-line arguments.  These should probably be stored in

	// e.g. Mandatory option:
	// numRepeats = convertTo<int>(args.Read("reps"));
	// An error will be thrown if there's no --reps=NNN option.

	// e.g. Optional option:
	// T1 = convertTo<double>(args.ReadWithDefault("T1","1.6"))

	// e.g. a basis function:
	// string basisFile = args.Read("basis"); // error if --basis=??? option is missing
	// basis = read_vest_fabber(basisFile); // error if file is missing

	// e.g. a boolean option (present or absent);
	// wantFriesWithThat = args.ReadBool("chips");


	// TODO: Read any commmand-line arguments and save them in the class's member variables.


	// Before you return, all of the implementation-specific variables in fwdmodel_custom.h
	// should have had values assigned to them.  If you didn't add any variables then you
	// don't have to do anything here!
	return;
}

void CustomFwdModel::Evaluate(const NEWMAT::ColumnVector& params, NEWMAT::ColumnVector& result) const
{
	assert(params.Nrows() == NumParams());

	// TODO:This needs to equal the number of timepoints in your data.
	// If it varies, you need to be able to calculate it from the command-line options.
	// This may mean adding a (redundant) --data-length=NNN option above.
	const int dataLength = 10; // Example value

	if (result.Nrows() != dataLength)
		result.ReSize(dataLength);

	// TODO: Your model here!

	// A very simple forward model: fitting a quadratic equation, with inputs at 1..10
	// Notice that NEWMAT vectors and matrices could from one, not zero!
	const double constantTerm = params(1);
	const double linearTerm = params(2);
	const double quadraticTerm = params(3);

	for (int i = 1; i <= dataLength; i++)
	{
		result(i) = constantTerm + i * linearTerm + i * i * quadraticTerm;
	}

	// Before you return, you should have assigned a predicted signal to "result".
	return;
}

int CustomFwdModel::NumParams() const
{
	// TODO: How many parameters does your model need?

	return 3; // Example function
}

void CustomFwdModel::NameParams(vector<string>& names) const
{
	// TODO: Name your model parameters here.  Should match NumParams() above!

	// Example function:
	names.push_back("Constant term");
	names.push_back("Linear term");
	names.push_back("Quadratic term");
}

string CustomFwdModel::ModelVersion() const
{
	return "$Id: fwdmodel_custom.cc,v 1.2 2012/01/13 12:00:59 adriang Exp $";
}

void CustomFwdModel::HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const
{
	Tracer_Plus tr("CustomFwdModel::HardcodedInitialDists");
	// Pick a safe starting point for your model fits, if no other input is provided.
	// The default one just uses a N(0,1e12) prior, and starts with all parameters set to zeroes:

	assert(prior.means.Nrows() == NumParams());

	prior.means = 0;
	prior.SetPrecisions(NEWMAT::IdentityMatrix(NumParams()) * 1e-12);
	posterior = prior;

	// Ask Adrian for help if you want to modify this!
	// You can override prior at the command line with the --fwd-initial-prior=VESTFILE option,
	// and you can set initialization points using the --continue-from-mvn=MVNFILE option.
}

void CustomFwdModel::Usage(std::ostream &stream)
{
	stream << "No model usage info available for --model=custom, yet." << endl;
}
