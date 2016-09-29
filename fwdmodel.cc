/*  fwdmodel.cc - base class for generic forward models

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"

#include "easylog.h"

#include <sstream> 
#include <memory>

FwdModel* FwdModel::NewFromName(const string& name)
{
	FwdModelFactory* factory = FwdModelFactory::GetInstance();
	FwdModel* model = factory->Create(name);
	if (model == NULL)
	{
		throw Invalid_option("Unrecognized forward model --model: " + name);
	}
	return model;
}

void FwdModel::UsageFromName(const string& name, std::ostream &stream)
{
	stream << "Usage information for model: " << name << endl << endl;
	std::auto_ptr<FwdModel> model(NewFromName(name));
	stream << model->GetDescription() << endl << endl << "Options: " << endl << endl;
	vector<OptionSpec> options;
	model->GetOptions(options);
	if (options.size() > 0)
	{
		for (vector<OptionSpec>::iterator iter = options.begin(); iter != options.end(); iter++)
		{
			stream << "  " << iter->name << " : " << iter->description << " "
					<< (iter->optional ? "(optional, default=" + iter->def : "(mandatory)") << endl;
		}
	}
	else {
		model->Usage(stream);
	}
}

string FwdModel::ModelVersion() const
{
	// You should overload this function in your FwdModel class

	// Something like this:
	// return " $ I d $ "; // without the spaces
	// CVS will automatically replace this with version information that looks
	// like this: $Id: fwdmodel.cc,v 1.38 2014/10/24 15:25:27 chappell Exp $
	// It should probably go in your .cc file, not the header.
	return "No version info available.";
}

void FwdModel::Usage(std::ostream &stream) const
{
	stream << "No usage information available" << endl;
}

bool FwdModel::Gradient(const ColumnVector& params, Matrix& grad) const
{
	// By default return false -> no gradient is supplied by this model
	return false;
}

void FwdModel::DumpParameters(const ColumnVector& params, const string& indent) const
{
	LOG << indent << "Parameters:" << endl;
	vector<string> names;
	NameParams(names);
	assert(names.size() == params.Nrows());

	for (int i = 1; i <= names.size(); i++)
		LOG << indent << "  " << names[i - 1] << " = " << params(i) << endl;

	LOG << indent << "Total of " << names.size() << " parameters" << endl;
}

void FwdModel::pass_in_coords(const ColumnVector& coords)
{
	coord_x = coords(1);
	coord_y = coords(2);
	coord_z = coords(3);
}

