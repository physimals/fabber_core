/*  fwdmodel.cc - base class for generic forward models

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

    Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"

#include <sstream> 
#include "easylog.h"
 
string FwdModel::ModelVersion() const
{
  return "No version info available.";
  // You should overload this function in your FwdModel class
  
  // Something like this:
  // return " $ I d $ "; // without the spaces
  // CVS will automatically replace this with version information that looks
  // like this: $Id: fwdmodel.cc,v 1.38 2014/10/24 15:25:27 chappell Exp $
  // It should probably go in your .cc file, not the header.
}

vector<string> FwdModel::GetUsage() const { 
  vector<string> usage;
  usage.push_back("General: undefined");
  usage.push_back("Required parameters: undefined");
  usage.push_back("Optional parameters: undefined");
  return usage;
}

int FwdModel::Gradient(const ColumnVector& params, Matrix& grad) const
{
  // by default return false -> no gradient is supplied by this model
  return false;
}

int FwdModel::NumOutputs() const
{
    ColumnVector params, result;
    params.ReSize(NumParams()); params = 1; // probably safer than 0
    Evaluate(params, result);
    return result.Nrows();
}

void FwdModel::DumpParameters(const ColumnVector& params, 
                              const string& indent) const
{
    LOG << indent << "Parameters:" << endl;
    vector<string> names;
    NameParams(names);
    
    for (int i = 1; i <= NumParams(); i++)
        LOG << indent << "  " << names[i-1] << " == " << params(i) << endl;
        
    LOG << indent << "Total of " << NumParams() << " parameters." << endl;
}

void FwdModel::pass_in_coords( const ColumnVector& coords )
{
  coord_x = coords(1);
  coord_y = coords(2);
  coord_z = coords(3);
}


