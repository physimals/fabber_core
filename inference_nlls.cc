#ifndef NO_NLLS

/* inference_nlls.cc - Non-Linear Least Squares class declarations

 Adrian Groves Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford */

/*  CCOPYRIGHT  */

#include "inference_nlls.h"

#include "easylog.h"
#include "fwdmodel.h"
#include "priors.h"
#include "rundata.h"
#include "tools.h"
#include "version.h"

#include <newmat.h>

#include <string>
#include <vector>

using namespace std;
using namespace MISCMATHS;
using namespace NEWMAT;
using fabber::MaskRows;

static int NUM_OPTIONS = 1;
static OptionSpec OPTIONS[] = {
    { "vb-init", OPT_BOOL, "Whether NLLS is being run in isolation or as a pre-step for VB",
        OPT_NONREQ, "" },
    { "lm", OPT_BOOL, "Whether to use LM convergence (default is L)", OPT_NONREQ, "" },
};

void NLLSInferenceTechnique::GetOptions(vector<OptionSpec> &opts) const
{
    InferenceTechnique::GetOptions(opts);
    for (int i = 0; i < NUM_OPTIONS; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

std::string NLLSInferenceTechnique::GetDescription() const
{
    return "Non-linear least squares inference technique.";
}

string NLLSInferenceTechnique::GetVersion() const
{
    return fabber_version();
}

InferenceTechnique *NLLSInferenceTechnique::NewInstance()
{
    return new NLLSInferenceTechnique();
}

void NLLSInferenceTechnique::Initialize(FabberRunData &args)
{
    InferenceTechnique::Initialize(args);
    LOG << "NLLSInferenceTechnique::Initialising" << endl;

    // Determine whether we use L (default) or LM convergence
    m_lm = args.GetBool("lm");
    LOG << "NLLSInferenceTechnique::Done initialising" << endl;
}

void NLLSInferenceTechnique::DoCalculations(FabberRunData &allData)
{
    m_ctxs.push_back(new ThreadContext(allData));
    ThreadContext *ctx = m_ctxs[0];

    // Check how many samples in time series (ignoring any masked time points)
    int Nsamples = ctx->m_origdata.Nrows() - ctx->m_masked_tpoints.size();

    // Loop over voxels. The result for each voxel is
    // stored as a MVN distribution for its parameters
    // in resultMVNs.
    for (unsigned int voxel = 1; voxel <= ctx->nvoxels; voxel++)
    {
        ColumnVector y = ctx->m_origdata.Column(voxel);
        ColumnVector vcoords = ctx->m_coords.Column(voxel);

        // Some models might want more information about the data
        ctx->PassModelData(voxel);
        IdentityMatrix I(ctx->m_num_params);

        // Create a cost function evaluator which will
        // measure the difference between the model
        // and the data
        NLLSCF costfn(y, ctx->m_model, ctx->m_masked_tpoints);

        // Set the convergence method
        // either Levenberg (L) or Levenberg-Marquardt (LM)
        NonlinParam nlinpar(ctx->m_num_params, NL_LM);
        if (!m_lm)
        {
            nlinpar.SetGaussNewtonType(LM_L);
        }

        // set ics from 'posterior'
        nlinpar.SetStartingEstimate(ctx->fwd_post[voxel-1].means);
        nlinpar.LogPar(true);
        nlinpar.LogCF(true);

        try
        {
            // Run the nonlinear optimizer
            // output variable is unused - unsure if nonlin has any effect
            nonlin(nlinpar, costfn);

#if 0
			LOG << "NLLSInferenceTechnique::The solution is: " << nlinpar.Par() << endl;
			LOG << "NLLSInferenceTechnique::and this is the process " << endl;
			for (int i=0; i<nlinpar.CFHistory().size(); i++)
			{
				LOG << " cf: " << (nlinpar.CFHistory())[i] <<endl;
			}
			for (int i=0; i<nlinpar.ParHistory().size(); i++)
			{
				LOG << (nlinpar.ParHistory())[i] << ": :";
			}
#endif
            // Get the new parameters
            ctx->fwd_post[voxel-1].means = nlinpar.Par();

            // Recenter linearized model on new parameters
            ctx->m_lin_model[voxel-1].ReCentre(ctx->fwd_post[voxel-1].means);
            Matrix J = ctx->m_lin_model[voxel-1].Jacobian();
            MaskRows(J, ctx->m_masked_tpoints);

            // Calculate the NLLS precision
            // This is (J'*J)/mse
            // The covariance is the inverse
            SymmetricMatrix nllsprec;
            double sqerr = costfn.cf(ctx->fwd_post[voxel-1].means);
            double mse = sqerr / (Nsamples - ctx->m_num_params);
            nllsprec << J.t() * J / mse;

            // Look for zero diagonal elements (implies parameter is not observable)
            // and set precision small, but non-zero - so that covariance can be calculated
            for (int i = 1; i <= nllsprec.Nrows(); i++)
            {
                if (nllsprec(i, i) < 1e-6)
                {
                    nllsprec(i, i) = 1e-6;
                }
            }
            ctx->fwd_post[voxel-1].SetPrecisions(nllsprec);
            ctx->fwd_post[voxel-1].GetCovariance();
        }
        catch (Exception &e)
        {
            LOG << "NLLSInferenceTechnique::NEWMAT Exception in this voxel:\n" << e.what() << endl;

            if (m_halt_bad_voxel)
                throw;

            LOG << "NLLSInferenceTechnique::Estimates in this voxel may be unreliable" << endl
                << "   (precision matrix will be set manually)" << endl
                << "   Going on to the next voxel" << endl;

            // output the results where we are
            ctx->fwd_post[voxel-1].means = nlinpar.Par();

            // recenter linearized model on new parameters
            ctx->m_lin_model[voxel-1].ReCentre(ctx->fwd_post[voxel-1].means);

            // precision matrix is probably singular so set manually
            ctx->fwd_post[voxel-1].SetPrecisions(I * 1e-12);
        }

        ctx->resultMVNs.at(voxel-1) = new MVNDist(ctx->fwd_post[voxel-1]);
    }
}

NLLSCF::NLLSCF(
    const NEWMAT::ColumnVector &pdata, const FwdModel *pm, std::vector<int> masked_tpoints)
    : m_data(MaskRows(pdata, masked_tpoints))
    , m_model(pm)
    , m_linear(pm)
    , m_masked_tpoints(masked_tpoints)
{
}

double NLLSCF::cf(const ColumnVector &p) const
{
    // p = parameters
    // data_pred = data predicted by model
    ColumnVector data_pred;
    m_model->EvaluateFabber(p, data_pred, "");
    data_pred = MaskRows(data_pred, m_masked_tpoints);

    // m_data = actual data. Find sum of squares of differences
    // between this and the model data using a scalar product.
    double cfv = ((m_data - data_pred).t() * (m_data - data_pred)).AsScalar();
    return cfv;
}

ReturnMatrix NLLSCF::grad(const ColumnVector &p) const
{
    // Create an initial zero gradient vector
    ColumnVector gradv(p.Nrows());

    // Need to recenter the linearised model to the current parameter values
    m_linear.ReCentre(p);
    Matrix J = m_linear.Jacobian();
    J = MaskRows(J, m_masked_tpoints);

    // Evaluate the model given the parameters
    ColumnVector data_pred;
    m_model->EvaluateFabber(p, data_pred, "");
    data_pred = MaskRows(data_pred, m_masked_tpoints);

    gradv = -2 * J.t() * (m_data - data_pred);
    gradv.Release();
    return (gradv);
}

boost::shared_ptr<BFMatrix> NLLSCF::hess(
    const ColumnVector &p, boost::shared_ptr<BFMatrix> iptr) const
{
    boost::shared_ptr<BFMatrix> hessm;

    if (iptr && iptr->Nrows() == (unsigned)p.Nrows() && iptr->Ncols() == (unsigned)p.Nrows())
    {
        hessm = iptr;
    }
    else
    {
        hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(), p.Nrows()));
    }

    // need to recenter the linearised model to the current parameter values
    m_linear.ReCentre(p);
    Matrix J = m_linear.Jacobian();
    J = MaskRows(J, m_masked_tpoints);
    Matrix hesstemp = 2 * J.t() * J; // Make the G-N approximation to the hessian

    //(*hessm) = J.t()*J;

    for (int i = 1; i <= p.Nrows(); i++)
    {
        for (int j = 1; j <= p.Nrows(); j++)
            hessm->Set(i, j, hesstemp(i, j));
    }

    return (hessm);
}

#endif
