/*  fwdmodel.cc - base class for generic forward models

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"

#include <sstream> 
#include <memory>
#include "easylog.h"

FwdModel* FwdModel::NewFromName(const string& name)
{
    FwdModelFactory* factory = FwdModelFactory::GetInstance();
    FwdModel* model = factory->Create(name);
    if (model == NULL) {
        throw Invalid_option("Unrecognized forward model --model: " + name);
    }
    return model;
}

void FwdModel::UsageFromName(const string& name, std::ostream &stream)
{
    std::auto_ptr<FwdModel> model(NewFromName(name));
    model->Usage(stream);
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

int FwdModel::Gradient(const ColumnVector& params, Matrix& grad) const
{
    // By default return false -> no gradient is supplied by this model
    return false;
}

int FwdModel::NumOutputs() const
{
    ColumnVector params, result;
    params.ReSize(NumParams());
    params = 1; // probably safer than 0
    Evaluate(params, result);
    return result.Nrows();
}

void FwdModel::DumpParameters(const ColumnVector& params, const string& indent) const
{
    LOG << indent << "Parameters:" << endl;
    vector<string> names;
    NameParams(names);

    for (int i = 1; i <= NumParams(); i++)
        LOG << indent << "  " << names[i - 1] << " == " << params(i) << endl;

    LOG << indent << "Total of " << NumParams() << " parameters." << endl;
}

void FwdModel::pass_in_coords(const ColumnVector& coords)
{
    coord_x = coords(1);
    coord_y = coords(2);
    coord_z = coords(3);
}

