/*  fwdmodel_flobs.cc - Does FLOBS

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */


#include "fwdmodel_flobs.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

string FlobsFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_flobs.cc,v 1.5 2012/01/13 12:00:59 adriang Exp $";
}

void FlobsFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("FlobsFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());
    
    // Set priors
    prior.means = 0;
    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    prior.SetPrecisions(precisions);
    
    posterior = prior;
}

void FlobsFwdModel::Initialize(ArgsType& args) 
{
    string basisFile = args.Read("basis"); 
    LOG_ERR( "    Reading basis functions: " << basisFile << endl );
    basis = read_vest(basisFile);
    LOG_ERR( "    Read " << basis.Ncols() << " basis functions of length " 
	     << basis.Nrows() << endl);

    basisFile = args.ReadWithDefault("nuisance","null"); 
    if (basisFile == "null")
      {
	nuisanceBasis.ReSize(basis.Nrows(), 0);
      }
    else if (basisFile == "offset")
      {
	nuisanceBasis.ReSize(basis.Nrows(), 1);
	nuisanceBasis = 1;
      }
    else
      {
	LOG_ERR( "    Reading nuisance basis functions: " 
		 << basisFile << endl );
	nuisanceBasis = read_vest(basisFile);
	LOG_ERR( "    Read " << nuisanceBasis.Ncols() 
		 << " nuisance basis functions of length "
		 << nuisanceBasis.Nrows() << endl);
      }

    if (nuisanceBasis.Nrows() != basis.Nrows())
      throw Invalid_option("Basis length mismatch!\n");
}

void FlobsFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    LOG << indent << "Scale = " << vec(1) << ", shape:\n1 (fixed)\n"
	<< vec.Rows(2, NumParams());
  // TODO: should dump nuisanceBasis too
}

void FlobsFwdModel::NameParams(vector<string>& names) const
{
    names.clear();
    
    for (int i = 1; i <= basis.Ncols(); i++)
      names.push_back("basis_" + stringify(i));

    for (int i = 1; i <= nuisanceBasis.Ncols(); i++)
      names.push_back("nuisance_"+stringify(i));
  
    assert(names.size() == (unsigned)NumParams()); 
}

//--- The Flobs 7 forward model - with polar co-ords --//


FactoryRegistration<FwdModelFactory,  Flobs7FwdModel> 
  Flobs7FwdModel::registration("flobs7");

void Flobs7FwdModel::Evaluate(const ColumnVector& params, 
                              ColumnVector& result) const {
  Tracer_Plus tr("Flobs7FwdModel::Evaluate");

  assert(params.Nrows() == NumParams());

  assert(basis.Ncols() == 2);
  ColumnVector beta = params.Rows(1,basis.Ncols());
  double betaBar = params(1);
  beta(1) = cos(beta(2));
  beta(2) = sin(beta(2));
  result = basis * beta * betaBar;
  // No nuisance stuff yet.
  assert(params.Nrows() == 2);

  return;
}

vector<string> Flobs7FwdModel::GetUsage() const { 
  vector<string> usage;
  usage.push_back("Required parameters:");
  usage.push_back("--basis=<basis_functions>");
  usage.push_back("--nuiscance=null|offset");
  usage.push_back("A scale parameter will automatically be added.");
  usage.push_back("Always use with the --flobs-prior-adjust option!!!");
  return usage;
}

FwdModel* Flobs7FwdModel::NewInstance() {
  return new Flobs7FwdModel();
}

//--- The Flobs 5 forward model --//

FactoryRegistration<FwdModelFactory,  Flobs5FwdModel> 
  Flobs5FwdModel::registration("flobs5");

void Flobs5FwdModel::Evaluate(const ColumnVector& params, 
                              ColumnVector& result) const {

  Tracer_Plus tr("Flobs5FwdModel::Evaluate");

  assert(params.Nrows() == NumParams());
  
  ColumnVector beta = params.Rows(1,basis.Ncols());
  double betaBar = params(1);
  beta(1) = 1.0; // fixed shape=1, scale=betaBar
  result = basis * beta * betaBar;
    
  beta = params.Rows(basis.Ncols()+1,params.Nrows());
  result += nuisanceBasis * beta;

  return; // answer is in the "result" vector
}

vector<string> Flobs5FwdModel::GetUsage() const { 
  vector<string> usage;
  usage.push_back("Required parameters:");
  usage.push_back("--basis=<basis_functions>");
  usage.push_back("--nuiscance=null|offset");
  usage.push_back("The first basis function will serve as the scaling factor (fixed shape==1)");
  return usage;
}

FwdModel* Flobs5FwdModel::NewInstance() {
  return new Flobs5FwdModel();
}
