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
  return "$Id: fwdmodel_flobs.cc,v 1.4 2008/04/03 13:22:39 adriang Exp $";
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
    //    if (useSeparateScale)
    //      {
    //	// Set informative initial posterior
    //	// Shape = first basis function, magnitude 1.
    //	posterior.means(1) = 1;
    //	posterior.means(NumParams()) = 1;
    //      }
}

void FlobsFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("FlobsFwdModel::Evaluate");

  assert(params.Nrows() == NumParams());
  
  //  if (useSeparateScale)
  //    {
  //      ColumnVector beta = params.Rows(1,basis.Ncols());
  //      //double betaBar = params(NumParams());
  //      result = basis * beta; // * betaBar; Now it's purely linear and ignore betaBar!
  //      assert(false); // TODO: update for nuisance parameters
  //    }

  if (usePolarCoordinates)
    {
      assert(basis.Ncols() == 2);
      ColumnVector beta = params.Rows(1,basis.Ncols());
      double betaBar = params(1);
      beta(1) = cos(beta(2));
      beta(2) = sin(beta(2));
      result = basis * beta * betaBar;

      // No nuisance stuff yet.
      assert(params.Nrows() == 2);
    }
  else
    { 
      ColumnVector beta = params.Rows(1,basis.Ncols());
      double betaBar = params(1);
      beta(1) = 1.0; // fixed shape=1, scale=betaBar
      result = basis * beta * betaBar;
      
      beta = params.Rows(basis.Ncols()+1,params.Nrows());
      result += nuisanceBasis * beta;
    }

  return; // answer is in the "result" vector
}

void FlobsFwdModel::ModelUsage()
{
  //  if (useSeparateScale)
    {
      cout << "Usage for --model=flobs6:\n"
	   << "  --basis=<basis_functions>\n"
	   << "  A scale parameter will automatically be added.\n"
	   << "  Always use with the --flobs-prior-adjust option!!! \n"
	   << "  (Future work: specify several scaling factors, for multi-event stimuli)\n\n";
    } 
    //  else 
    {
      cout << "Usage for --model=flobs5:\n"
	   << "  --basis=<basis_functions>\n"
	   << "  The first basis function will serve as the scaling factor (fixed shape==1)\n";
    }
    // TODO: add --nuisance= option
}

FlobsFwdModel::FlobsFwdModel(ArgsType& args, bool polar) 
//  : useSeparateScale(sepScale)
  : usePolarCoordinates(polar)
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
  //  if (useSeparateScale)
  //    LOG << indent << "Scale = " << vec(NumParams()) << ", shape:\n"
  //	<< vec.Rows(1, NumParams()-1);
  //  else
    LOG << indent << "Scale = " << vec(1) << ", shape:\n1 (fixed)\n"
	<< vec.Rows(2, NumParams());
  // TODO: should dump nuisanceBasis too
}

void FlobsFwdModel::NameParams(vector<string>& names) const
{
    names.clear();
    
    for (int i = 1; i <= basis.Ncols(); i++)
      names.push_back("basis_" + stringify(i));
    
    //    if (useSeparateScale)
    //      names.push_back("scale");

    for (int i = 1; i <= nuisanceBasis.Ncols(); i++)
      names.push_back("nuisance_"+stringify(i));
  
    assert(names.size() == (unsigned)NumParams()); 
}
