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
static OptionSpec OPTIONS[] = { { "basis", OPT_MATRIX, "Design matrix", OPT_REQ, "" } };

FwdModel *LinearFwdModel::NewInstance() { return new LinearFwdModel(); }
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

string LinearFwdModel::ModelVersion() const { return fabber_version(); }
void LinearFwdModel::Initialize(FabberRunData &args)
{
    FwdModel::Initialize(args);
    string designFile = args.GetString("basis");
    LOG << "LinearFwdModel::Reading design file: " << designFile << endl;
    m_jacobian = read_matrix_file(designFile);

    const int Ntimes = m_jacobian.Nrows();
    const int Nbasis = m_jacobian.Ncols();

    LOG << "LinearFwdModel::Loaded " << m_jacobian.Ncols() << " basis functions of length "
        << Ntimes << endl;

    m_centre.ReSize(Nbasis);
    m_centre = 0;
    m_offset.ReSize(Ntimes);
    m_offset = 0;

    if (args.GetBool("add-ones-regressor"))
    {
        // Add an additional 'parameter' whose timeseries if constant
        LOG << "LinearFwdModel::Plus an additional regressor of all ones\n";
        ColumnVector ones(Ntimes);
        ones = 1.0;
        m_jacobian = m_jacobian | ones;
        m_centre.ReSize(Nbasis + 1);
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

void LinearFwdModel::EvaluateModel(
    const ColumnVector &params, ColumnVector &result, const std::string &key) const
{
    result = m_jacobian * (params - m_centre) + m_offset;
}

int LinearFwdModel::NumParams() const { return m_centre.Nrows(); }
void LinearFwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    // Completely generic names.
    for (int i = 1; i <= NumParams(); i++)
        names.push_back("Parameter_" + stringify(i));
}

LinearizedFwdModel::LinearizedFwdModel(const FwdModel *model)
    : m_model(model)
{
    SetLogger(model->GetLogger());
}

LinearizedFwdModel::LinearizedFwdModel(const LinearizedFwdModel &from)
    : LinearFwdModel(from)
    , m_model(from.m_model)
{
    SetLogger(from.GetLogger());
}

void LinearizedFwdModel::HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const
{
    assert(m_model);
    m_model->HardcodedInitialDists(prior, posterior);
}

void LinearizedFwdModel::NameParams(std::vector<std::string> &names) const
{
    assert(m_model);
    m_model->NameParams(names);
}

std::string LinearizedFwdModel::ModelVersion()
{
    assert(m_model != NULL);
    return m_model->ModelVersion();
}

void LinearizedFwdModel::ReCentre(const ColumnVector &about)
{
    assert(about == about); // isfinite

    // Store new centre & offset
    m_centre = about;

    m_model->EvaluateFabber(m_centre, m_offset);
    if (0 * m_offset != 0 * m_offset)
    {
        LOG_ERR("LinearizedFwdModel::about:\n" << about);
        LOG_ERR("LinearizedFwdModel::m_offset:\n" << m_offset.t());
        throw FabberInternalError(
            "LinearizedFwdModel::ReCentre: Non-finite values found in offset");
    }

    // Calculate the Jacobian numerically.  jacobian is len(y)-by-len(m)
    m_jacobian.ReSize(m_offset.Nrows(), m_centre.Nrows());
    // m_jacobian = 0.0/0.0; // fill with NaNs to check

    // Try and get the gradient matrix (Jacobian) from the model first
    // FIXME this is broken when transforms are used
    // int gradfrommodel = m_model->Gradient(m_centre, m_jacobian);

    // If the gradient is not supported by the model, use
    // numerical differentiation to calculate it.
    if (true)
    {
        ColumnVector centre2, centre3;
        ColumnVector offset2, offset3;
        for (int i = 1; i <= m_centre.Nrows(); i++)
        {
            double delta = m_centre(i) * 1e-5;
            if (delta < 0)
                delta = -delta;
            if (delta < 1e-10)
                delta = 1e-10;

            // Take derivative numerically
            centre3 = m_centre;
            centre2 = m_centre;
            centre2(i) += delta;
            centre3(i) -= delta;
            m_model->EvaluateFabber(centre2, offset2);
            m_model->EvaluateFabber(centre3, offset3);
            m_jacobian.Column(i) = (offset2 - offset3) / (centre2(i) - centre3(i));
        }
    }

    if (0 * m_jacobian != 0 * m_jacobian)
    //    if(true)
    {
        LOG << "LinearizedFwdModel::jacobian:\n" << m_jacobian;
        LOG << "LinearizedFwdModel::about':\n" << about.t();
        LOG << "LinearizedFwdModel::offset':\n" << m_offset.t();
        throw FabberInternalError(
            "LinearizedFwdModel::ReCentre: Non-finite values found in jacobian");
    }
}

void LinearizedFwdModel::DumpParameters(const ColumnVector &vec, const string &indent) const
{
    //    LOG << indent << "This is what the nonlinear model has to say:" << endl;
    m_model->DumpParameters(vec, indent);
}
