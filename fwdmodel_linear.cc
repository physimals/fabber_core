/*  fwdmodel_linear.cc - Linear forward model and related classes

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_linear.h"

#include "dist_mvn.h"
#include "easylog.h"
#include "rundata.h"
#include "tools.h"
#include "version.h"

#include <newmatio.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using fabber::read_matrix_file;
using namespace std;
using namespace NEWMAT;

static int NUM_OPTIONS = 1;
static OptionSpec OPTIONS[] = {
    { "basis", OPT_MATRIX, "File containing design matrix", OPT_REQ, "" }
};

FwdModel *LinearFwdModel::NewInstance()
{
    return new LinearFwdModel();
}

void LinearFwdModel::GetOptions(std::vector<OptionSpec> &opts) const
{
    for (int i = 0; i < NUM_OPTIONS; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

std::string LinearFwdModel::GetDescription() const
{
    return "Model in which output is a linear combination of input parameters";
}

string LinearFwdModel::ModelVersion() const
{
    return fabber_release_version();
}

void LinearFwdModel::Initialize(FabberRunData &args)
{
    FwdModel::Initialize(args);
    string designFile = args.GetString("basis");
    LOG << "LinearFwdModel::Reading design file: " << designFile << endl;
    jacobian = read_matrix_file(designFile);

    const int Ntimes = jacobian.Nrows();
    const int Nbasis = jacobian.Ncols();

    LOG << "LinearFwdModel::Loaded " << jacobian.Ncols() << " basis functions of length " << Ntimes << endl;

    centre.ReSize(Nbasis);
    centre = 0;
    offset.ReSize(Ntimes);
    offset = 0;

    if (args.GetBool("add-ones-regressor"))
    {
        // Add an additional 'parameter' whose timeseries if constant
        LOG << "LinearFwdModel::Plus an additional regressor of all ones\n";
        ColumnVector ones(Ntimes);
        ones = 1.0;
        jacobian = jacobian | ones;
        centre.ReSize(Nbasis + 1);
    }

    //
    // Warning: Nbasis is now wrong!
}

void LinearFwdModel::HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const
{
    assert(prior.means.Nrows() == NumParams());

    // Set default for each parameter to a mean of 0 and close-to-zero precision.
    prior.means = 0;
    prior.SetPrecisions(IdentityMatrix(NumParams()) * 1e-12);
    posterior = prior;
}

void LinearFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    result = jacobian * (params - centre) + offset;
}

int LinearFwdModel::NumParams() const
{
    return centre.Nrows();
}

void LinearFwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    // Completely generic names.
    for (int i = 1; i <= NumParams(); i++)
        names.push_back("Parameter_" + stringify(i));
}

void LinearizedFwdModel::ReCentre(const ColumnVector &about)
{
    assert(about == about); // isfinite

    // Store new centre & offset
    centre = about;
    fcn->Evaluate(centre, offset);
    if (0 * offset != 0 * offset)
    {
        LOG_ERR("LinearizedFwdModel::about:\n"
            << about);
        LOG_ERR("LinearizedFwdModel::offset:\n"
            << offset.t());
        throw FabberInternalError("LinearizedFwdModel::ReCentre: Non-finite values found in offset");
    }

    // Calculate the Jacobian numerically.  jacobian is len(y)-by-len(m)
    jacobian.ReSize(offset.Nrows(), centre.Nrows());
    // jacobian = 0.0/0.0; // fill with NaNs to check

    // Try and get the gradient matrix (Jacobian) from the model first
    int gradfrommodel = fcn->Gradient(centre, jacobian);

    // If the gradient is not supported by the model, use
    // numerical differentiation to calculate it.
    if (!gradfrommodel)
    {
        ColumnVector centre2, centre3;
        ColumnVector offset2, offset3;
        for (int i = 1; i <= centre.Nrows(); i++)
        {
            double delta = centre(i) * 1e-5;
            if (delta < 0)
                delta = -delta;
            if (delta < 1e-10)
                delta = 1e-10;

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

    if (0 * jacobian != 0 * jacobian)
    {
        LOG << "LinearizedFwdModel::jacobian:\n"
            << jacobian;
        LOG << "LinearizedFwdModel::about':\n"
            << about.t();
        LOG << "LinearizedFwdModel::offset':\n"
            << offset.t();
        throw FabberInternalError("LinearizedFwdModel::ReCentre: Non-finite values found in jacobian");
    }
}

void LinearizedFwdModel::DumpParameters(const ColumnVector &vec, const string &indent) const
{
    //    LOG << indent << "This is what the nonlinear model has to say:" << endl;
    fcn->DumpParameters(vec, indent);
}
