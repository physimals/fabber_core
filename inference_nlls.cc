/* inference_nlls.h - Non-Linear Least Squares class declarations

 Adrian Groves Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford */

/*  CCOPYRIGHT  */

#include "dataset.h"
#include "inference_nlls.h"

InferenceTechnique* NLLSInferenceTechnique::NewInstance()
{
    return new NLLSInferenceTechnique();
}

void NLLSInferenceTechnique::Initialize(FwdModel* fwd_model, FabberRunData& args)
{
    Tracer_Plus tr("NLLSInferenceTechnique::Initialize");

    LOG << "Initialising NLLS method" << endl;
    model = fwd_model;

    // Determine whether NLLS is being run in isolation or as a pre-step for VB
    // This alters what we do if result is ill conditioned
    vbinit = args.GetBool("vb-init");

    // Initialize the model with MVN distributions for its parameters
    MVNDist* loadPosterior = new MVNDist(model->NumParams());
    MVNDist* junk = new MVNDist(model->NumParams());
    model->HardcodedInitialDists(*junk, *loadPosterior);

    // Option to load a 'posterior' which will allow the setting of intial parameter estimates for NLLS
    string filePosterior = args.GetStringDefault("fwd-inital-posterior", "modeldefault");
    if (filePosterior != "modeldefault") {
        LOG << "File posterior" << endl;
        loadPosterior->Load(filePosterior);
    }

    initialFwdPosterior = loadPosterior;
    loadPosterior = NULL;

    // Determine whether we use L (default) or LM convergence
    lm = args.GetBool("lm");
    LOG << "Done initialising" << endl;
}

void NLLSInferenceTechnique::DoCalculations(FabberRunData& allData)
{
    Tracer_Plus tr("NLLSInferenceTechnique::DoCalculations");

    // Get basic voxel data
    const Matrix& data = allData.GetMainVoxelData();
    const Matrix & coords = allData.GetVoxelCoords();
    unsigned int Nvoxels = data.Ncols();

    // pass in some (dummy) data/coords here just in case the model relies upon it
    // use the first voxel values as our dummies

    model->pass_in_data(data.Column(1));
    model->pass_in_coords(coords.Column(1));

    // Check how many samples in time series - should
    // be same as model outputs
    int Nsamples = data.Nrows();
    if (data.Nrows() != model->NumOutputs()) {
        throw Invalid_option("Data length (" + stringify(data.Nrows()) + ") does not match model's output length ("
                + stringify(model->NumOutputs()) + ")!");
    }

    // Loop over voxels. The result for each voxel is
    // stored as a MVN distribution for its parameters
    // in resultMVNs.
    for (unsigned int voxel = 1; voxel <= Nvoxels; voxel++) {
        ColumnVector y = data.Column(voxel);
        ColumnVector vcoords = coords.Column(voxel);

        // Some models might want more information about the data
        model->pass_in_data(y);
        model->pass_in_coords(vcoords);

        LinearizedFwdModel linear(model);

        // FIXME should be a single sensible way to get the
        // number of model parameters!
        int Nparams = initialFwdPosterior->GetSize();

        // FIXME how about a ctor for MVNDist which takes a size?
        MVNDist fwdPosterior;
        fwdPosterior.SetSize(Nparams);

        IdentityMatrix I(Nparams);

        // Create a cost function evaluator which will
        // measure the difference between the model
        // and the data
        NLLSCF costfn(y, model);

        // Set the convergence method
        // either Levenberg (L) or Levenberg-Marquardt (LM)
        NonlinParam nlinpar(Nparams, NL_LM);
        if (!lm) {
            nlinpar.SetGaussNewtonType(LM_L);
        }

        // set ics from 'posterior'
        ColumnVector nlinics = initialFwdPosterior->means;
        nlinpar.SetStartingEstimate(nlinics);
        nlinpar.LogPar(true);
        nlinpar.LogCF(true);

        try {
            // Run the nonlinear optimizer
            // output variable status is unused - unsure if nonlin has any effect
            NonlinOut status = nonlin(nlinpar, costfn);

#if 0
            LOG << "The solution is: " << nlinpar.Par() << endl;
            LOG << "and this is the process " << endl;
            for (int i=0; i<nlinpar.CFHistory().size(); i++) {
                LOG << " cf: " << (nlinpar.CFHistory())[i] <<endl;
            }
            for (int i=0; i<nlinpar.ParHistory().size(); i++) {
                LOG << (nlinpar.ParHistory())[i] << ": :";
            }
#endif
            // Get the new parameters
            fwdPosterior.means = nlinpar.Par();

            // Recenter linearized model on new parameters
            linear.ReCentre(fwdPosterior.means);
            const Matrix& J = linear.Jacobian();

            // Calculate the NLLS precision
            // This is (J'*J)/mse
            // The covariance is the inverse
            SymmetricMatrix nllsprec;
            double sqerr = costfn.cf(fwdPosterior.means);
            double mse = sqerr / (Nsamples - Nparams);
            nllsprec << J.t() * J / mse;

            // Look for zero diagonal elements (implies parameter is not observable)
            // and set precision small, but non-zero - so that covariance can be calculated
            for (int i = 1; i <= nllsprec.Nrows(); i++) {
                if (nllsprec(i, i) < 1e-6) {
                    nllsprec(i, i) = 1e-6;
                }
            }
            fwdPosterior.SetPrecisions(nllsprec);
            fwdPosterior.GetCovariance();
        } catch (Exception &e) {
            LOG << "   NEWMAT Exception in this voxel:\n" << e.what() << endl;

            //if (haltOnBadVoxel) throw;

            LOG << "   Estimates in this voxel may be unreliable" << endl
                    << "   (precision matrix will be set manually)" << endl << "   Going on to the next voxel" << endl;

            // output the results where we are
            fwdPosterior.means = nlinpar.Par();

            // recenter linearized model on new parameters
            linear.ReCentre(fwdPosterior.means);

            // precision matrix is probably singular so set manually
            fwdPosterior.SetPrecisions(I * 1e-12);
        }

        resultMVNs.push_back(new MVNDist(fwdPosterior));
        assert(resultMVNs.size() == voxel);
    }
}

NLLSInferenceTechnique::~NLLSInferenceTechnique()
{
}

double NLLSCF::cf(const ColumnVector& p) const
{
    Tracer_Plus tr("NLLSCF::cf");

    // p = parameters
    // yhat = data predicted by model
    ColumnVector yhat;
    model->Evaluate(p, yhat);

    // y = actual data. Find sum of squares of differences
    // between this and the model data using a scalar product.
    double cfv = ((y - yhat).t() * (y - yhat)).AsScalar();
    return (cfv);
}

ReturnMatrix NLLSCF::grad(const ColumnVector& p) const
{
    Tracer_Plus tr("NLLSCF::grad");

    // Create an initial zero gradient vector
    ColumnVector gradv(p.Nrows());
    gradv = 0.0;

    // Need to recenter the linearised model to the current parameter values
    linear.ReCentre(p);
    const Matrix& J = linear.Jacobian();

#if 0 // FIXME Old code?
    //this is g(w) i.e. model evaluated at current parameters?
    //const ColumnVector gm = linear.Offset();
#endif
    // Evaluate the model given the parameters
    ColumnVector yhat;
    model->Evaluate(p, yhat);

    gradv = -2 * J.t() * (y - yhat);
    gradv.Release();
    return (gradv);
}

boost::shared_ptr<BFMatrix> NLLSCF::hess(const ColumnVector& p, boost::shared_ptr<BFMatrix> iptr) const
{
    Tracer_Plus tr("NLLSCF::hess");

    boost::shared_ptr<BFMatrix> hessm;

    if (iptr && iptr->Nrows() == (unsigned) p.Nrows() && iptr->Ncols() == (unsigned) p.Nrows()) {
        hessm = iptr;
    }
    else {
        hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(), p.Nrows()));
    }

    // need to recenter the linearised model to the current parameter values
    linear.ReCentre(p);
    const Matrix& J = linear.Jacobian();
    Matrix hesstemp = 2 * J.t() * J; //Make the G-N approximation to the hessian

    //(*hessm) = J.t()*J;

    for (int i = 1; i <= p.Nrows(); i++) {
        for (int j = 1; j <= p.Nrows(); j++)
            hessm->Set(i, j, hesstemp(i, j));
    }

    return (hessm);
}

