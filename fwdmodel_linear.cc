/*  fwdmodel_linear.cc - Linear forward model and related classes

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_linear.h"
#include <iostream>
#include "newmatio.h"
#include <stdexcept>
#include "easylog.h"
#include "newimage/newimageall.h"
using namespace NEWIMAGE;

string LinearFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_linear.cc,v 1.18 2011/08/04 13:36:35 chappell Exp $";
}

void LinearFwdModel::ModelUsage()
{
  cout << "\nUsage info for --model=linear:\n"
       << "Required options:\n"
       << "--basis=<design_file>\n"
    ;
}

LinearFwdModel::LinearFwdModel(ArgsType& args)
{
  Tracer_Plus tr("LinearFwdModel::LinearFwdModel(args)");
  string designFile = args.Read("basis");
  LOG_ERR("    Reading design file: " << designFile << endl);
  jacobian = read_vest(designFile);

  const int Ntimes = jacobian.Nrows();
  const int Nbasis = jacobian.Ncols();

  LOG_ERR("      Loaded " << jacobian.Ncols() 
	  << " basis functions of length " << Ntimes << endl);

  centre.ReSize(Nbasis); 
  centre = 0;
  offset.ReSize(Ntimes);
  offset = 0;

  if (args.ReadBool("add-ones-regressor"))
    {
      LOG_ERR("      Plus an additional regressor of all ones\n");
      ColumnVector ones(Ntimes); ones = 1.0;
      jacobian = jacobian | ones;
      centre.ReSize(Nbasis+1);
    }
  // Warning: Nbasis is now wrong!
}

void LinearFwdModel::HardcodedInitialDists(MVNDist& prior, 
					   MVNDist& posterior) const
{
  Tracer_Plus tr("LinearFwdModel::HardcodedInitialDists");
  assert(prior.means.Nrows() == NumParams());

  prior.means = 0;
  prior.SetPrecisions(IdentityMatrix(NumParams()) * 1e-12);
  posterior = prior;

}

void LinearFwdModel::Evaluate(const ColumnVector& params,
	               		    ColumnVector& result) const
{
  result = jacobian * (params - centre) + offset;
}

void LinearizedFwdModel::ReCentre(const ColumnVector& about)
{
  Tracer_Plus tr("LinearizedFwdModel::ReCentre");
  assert(about == about); // isfinite

  // Store new centre & offset
  centre = about;
  fcn->Evaluate(centre, offset);
  if (0*offset != 0*offset) 
    {
      LOG_ERR("about:\n" << about);
      LOG_ERR("offset:\n" << offset.t());
      throw overflow_error("ReCentre: Non-finite values found in offset");
    }

  // Calculate the Jacobian numerically.  jacobian is len(y)-by-len(m)
  jacobian.ReSize(offset.Nrows(), centre.Nrows());
  // jacobian = 0.0/0.0; // fill with NaNs to check
  
  // try and get the gradient from the model first
  int gradfrommodel=false;
  gradfrommodel = fcn->Gradient(centre,jacobian);

  if (!gradfrommodel) {
  ColumnVector centre2, centre3;
  ColumnVector offset2, offset3;
  for (int i = 1; i <= centre.Nrows(); i++)
    {
      double delta = centre(i) * 1e-5;
      if (delta<0) delta = -delta;
      if (delta<1e-10) delta = 1e-10;

      // Take derivative numerically
      centre3 = centre;
      centre2 = centre;
      centre2(i) += delta;
      centre3(i) -= delta;
      fcn->Evaluate(centre2, offset2);
      fcn->Evaluate(centre3, offset3);
      jacobian.Column(i) = (offset2 - offset3) / (centre2(i) - centre3(i));

      /*
if (i==4)
{LOG << "centre2 -centre3== \n" << 1e10*(centre2-centre3) << endl;
LOG << "offset2-offset3 == \n" << offset2(33)-offset3(33) << endl;
LOG << "offset2-offset3 == \n" << float(offset2(33)-offset3(33)) << endl;
LOG << "offset2-offset3 == \n" << double(offset2(33)-offset3(33)) << endl;
LOG << "Jac 33,4 == " << jacobian(33,4) << endl;
}
//*/
    }   
  }

  if (0*jacobian != 0*jacobian) 
    {
      LOG << "jacobian:\n" << jacobian;
      LOG << "about':\n" << about.t();
      LOG << "offset':\n" << offset.t();    
      throw overflow_error("ReCentre: Non-finite values found in jacobian");
    }
}

void LinearFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    LOG << indent << "Parameters mean nothing to me!  "
                   << "I am a mere linear model." << endl;
    LOG << indent << "  Vector: " << vec.t() << endl;
}

void LinearFwdModel::NameParams(vector<string>& names) const
{
    names.clear();
 
    // Completely generic names.
    for (int i = 1; i <= NumParams(); i++)
        names.push_back("Parameter_" + stringify(i));
}

void LinearizedFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
//    LOG << indent << "This is what the nonlinear model has to say:" << endl;
    fcn->DumpParameters(vec, indent);
}
