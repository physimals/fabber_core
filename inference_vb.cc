/*  inference_vb.cc - Variational Bayes with optional spatial smoothing

 Adrian Groves and Matthew Webster, FMRIB Image Analysis Group

 Copyright (C) 2007-2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_vb.h"

#include "convergence.h"
#include "easylog.h"
#include "priors.h"
#include "run_context.h"
#include "tools.h"
#include "version.h"
#include "fabber_threads.h"

#include <miscmaths/miscmaths.h>
#include <newmatio.h>

#include <math.h>
#include <algorithm>

using MISCMATHS::sign;

static OptionSpec OPTIONS[] = {
    { "noise", OPT_STR, "Noise model to use (white or ar1)", OPT_REQ, "" },
    { "convergence", OPT_STR, "Name of method for detecting convergence", OPT_NONREQ, "maxits" },
    { "max-iterations", OPT_STR,
        "number of iterations of VB to use with the maxits convergence detector", OPT_NONREQ,
        "10" },
    { "min-fchange", OPT_STR,
        "When using the fchange convergence detector, the change in F to stop at", OPT_NONREQ,
        "10" },
    { "max-trials", OPT_STR, "When using the trial mode convergence detector, the maximum number "
                             "of trials after an initial reduction in F",
        OPT_NONREQ, "10" },
    { "print-free-energy", OPT_BOOL, "Output the free energy", OPT_NONREQ, "" },
    { "mcsteps", OPT_INT, "Number of motion correction steps", OPT_NONREQ, "0" },
    { "continue-from-mvn", OPT_MVN, "Continue previous run from output MVN files", OPT_NONREQ, "" },
    { "output-only", OPT_BOOL, "Skip model fitting, just output requested data based on supplied "
                               "MVN. Can only be used with continue-from-mvn",
        OPT_NONREQ, "" },
    { "noise-initial-prior", OPT_MATRIX, "MVN of initial noise prior", OPT_NONREQ, "" },
    { "noise-initial-posterior", OPT_MATRIX, "MVN of initial noise posterior", OPT_NONREQ, "" },
    { "noise-pattern", OPT_STR, "repeating pattern of noise variances for each point (e.g. 12 "
                                "gives odd and even data points different variances)",
        OPT_NONREQ, "1" },
    { "PSP_byname<n>", OPT_STR, "Name of model parameter to use image prior", OPT_NONREQ, "" },
    { "PSP_byname<n>_type", OPT_STR, "Type of image prior to use for parameter <n> - I=image prior",
        OPT_NONREQ, "" },
    { "PSP_byname<n>_image", OPT_IMAGE, "Image prior for parameter <n>", OPT_NONREQ, "" },
    { "PSP_byname<n>_prec", OPT_FLOAT, "Precision to apply to image prior for parameter <n>",
        OPT_NONREQ, "" },
    { "PSP_byname<n>_transform", OPT_STR, "Transform to apply to parameter <n>", OPT_NONREQ, "" },
    { "allow-bad-voxels", OPT_BOOL,
        "Continue if numerical error found in a voxel, rather than stopping", OPT_NONREQ, "" },
    { "ar1-cross-terms", OPT_STR, "For AR1 noise, type of cross-linking (dual, same or none)",
        OPT_NONREQ, "dual" },
    { "spatial-dims", OPT_INT, "Number of spatial dimensions", OPT_NONREQ, "3" },
    { "spatial-speed", OPT_STR, "Restrict speed of spatial smoothing", OPT_NONREQ, "-1" },
    { "distance-measure", OPT_STR, "", OPT_NONREQ, "dist1" },
    { "param-spatial-priors", OPT_STR,
        "Type of spatial priors for each parameter, as a sequence of characters. "
        "N=nonspatial, M=Markov random field, P=Penny, A=ARD",
        OPT_NONREQ, "N+" },
    { "update-spatial-prior-on-first-iteration", OPT_BOOL, "", OPT_NONREQ, "" },
    { "locked-linear-from-mvn", OPT_MVN, "MVN file containing fixed centres for linearization",
        OPT_NONREQ, "" },
    { "" },
};

void Vb::GetOptions(vector<OptionSpec> &opts) const
{
    InferenceTechnique::GetOptions(opts);

    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string Vb::GetDescription() const
{
    return "Variational Bayes inference technique";
}
string Vb::GetVersion() const
{
    return fabber_version();
}
InferenceTechnique *Vb::NewInstance()
{
    return new Vb();
}
void Vb::Initialize(FabberRunData &rundata)
{
    InferenceTechnique::Initialize(rundata);

    // Motion correction related setup - by default no motion correction
    m_num_mcsteps = convertTo<int>(rundata.GetStringDefault("mcsteps", "0"));

    int m_spatial_dims = rundata.GetIntDefault("spatial-dims", 3, 0, 3);
    if (m_spatial_dims == 1)
    {
        WARN_ONCE("spatial-dims=1 is weird... hope you're just testing!");
    }
    else if (m_spatial_dims == 2)
    {
        WARN_ONCE("spatial-dims=2 may not work the way you expect");
    }
}

bool Vb::IsSpatial(ThreadContext *ctx, FabberRunData &rundata) const
{
    if (rundata.GetString("method") == "spatialvb")
    {
        return true;
    }
    else if (ctx)
    {
        // Really clunky way to detect if any spatial priors have been specified
        vector<Parameter> params;
        ctx->m_model->GetParameters(rundata, params);
        for (vector<Parameter>::iterator iter = params.begin(); iter != params.end(); iter++)
        {
            switch (iter->prior_type)
            {
            case PRIOR_SPATIAL_M:
            case PRIOR_SPATIAL_m:
            case PRIOR_SPATIAL_P:
            case PRIOR_SPATIAL_p:
                return true;
            }
        }
    }
    return false;
}

void Vb::DoCalculations(FabberRunData &rundata)
{
    int nthreads = rundata.GetIntDefault("nthreads", 1);
    int nvox = rundata.GetMainVoxelData().Ncols();
    nthreads = std::min(nthreads, nvox);
    nthreads = std::max(nthreads, 1);

    if (rundata.GetBool("output-only"))
    {
        m_ctxs.push_back(new VbThreadContext(rundata));
        // Do no calculations - now we have set resultMVNs we can finish
        LOG << "Vb::DoCalculations output-only set - not performing any calculations" << endl;
    }
    else if (IsSpatial(NULL, rundata)) // FIXME
    {
        // Make the neighbours[] lists if required
        
        std::vector<std::vector<int> > neighbours;
        std::vector<std::vector<int> > neighbours2;

        //if (m_prior_types_str.find_first_of("mMpP") != string::npos)
        if (true) // FIXME
        {
            int spatial_dims = rundata.GetIntDefault("spatial-dims", 3, 0, 3);
            CalcNeighbours(rundata.GetVoxelCoords(), spatial_dims, neighbours, neighbours2);
        }
        m_ctxs.push_back(new SpatialVbThreadContext(rundata, neighbours, neighbours2));
        ThreadContext *ctx = m_ctxs[0];
        ctx->Start(ctx);
    }
    else
    {
        int chunk = nvox / nthreads;
        int first = 1;
        vector<THREAD_ID> ids;
        for (int i=0; i<nthreads; i++) {
            if (i == (nthreads - 1))
            {
                // The last thread picks up all the remaining voxels
                chunk = 1 + nvox - first;
            }
            VbThreadContext *ctx = new VbThreadContext(rundata, i, nthreads, first, chunk);
            m_ctxs.push_back(ctx);
            first += chunk;
            THREAD_ID id;
            int status = create_thread(&id, &VbThreadContext::Start, (void *)ctx);
            if (status != 0) {
                throw FabberInternalError("Failed to create thread");
            } 
            cerr << "Started worker " << (i+1) << "/" << nthreads << endl;
            ids.push_back(id);
        }
        for (int i=0; i<nthreads; i++) {
            join_thread(ids[i]);
        }
    }
}

void Vb::SaveResults(FabberRunData &rundata) const
{
    InferenceTechnique::SaveResults(rundata);

    LOG << "Vb::Preparing to save results..." << endl;
    int total_nvoxels = 0;
    for (int i=0; i<m_ctxs.size(); i++) {
        total_nvoxels += m_ctxs[i]->nvoxels;
    }

    if (rundata.GetBool("save-noise-mean") | rundata.GetBool("save-noise-std"))
    {
        if (m_ctxs[0]->m_noise_params > 0)
        {
            LOG << "Vb::Writing noise" << endl;
            Matrix noiseMean, noiseStd;
            noiseMean.ReSize(m_ctxs[0]->m_noise_params, total_nvoxels);
            noiseStd.ReSize(m_ctxs[0]->m_noise_params, total_nvoxels);

            for (int c=0; c<m_ctxs.size(); c++) {
                ThreadContext *ctx = m_ctxs[c];
                for (int vox = 1; vox <= ctx->nvoxels; vox++)
                {
                    for (int i = 1; i <= ctx->m_noise_params; i++)
                    {
                        noiseStd(i, ctx->start_voxel + vox - 1) = sqrt(
                            ctx->resultMVNs[vox - 1]->GetCovariance()(i + ctx->m_num_params, i + ctx->m_num_params));
                        noiseMean(i, ctx->start_voxel + vox - 1) = ctx->resultMVNs[vox - 1]->means(i + ctx->m_num_params);
                    }
                }
            }

            // FIXME was this being saved before? Should it be?
            if (rundata.GetBool("save-noise-mean"))
                rundata.SaveVoxelData("noise_means", noiseMean);
            if (rundata.GetBool("save-noise-std"))
                rundata.SaveVoxelData("noise_stdevs", noiseStd);
        }
    }

    // Save the Free Energy estimates
    if (rundata.GetBool("save-free-energy") && !m_ctxs[0]->resultFs.empty())
    {
        LOG << "Vb::Writing free energy" << endl;
        Matrix freeEnergy;
        freeEnergy.ReSize(1, total_nvoxels);
        
        for (int c=0; c<m_ctxs.size(); c++) {
            ThreadContext *ctx = m_ctxs[c];
            for (int vox = 1; vox <= ctx->nvoxels; vox++)
            {
                freeEnergy(1, ctx->start_voxel + vox - 1) = ctx->resultFs.at(vox - 1);
            }
        }
        rundata.SaveVoxelData("freeEnergy", freeEnergy);
    }
    else
    {
        LOG << "Vb::Free energy wasn't recorded, so no freeEnergy data saved" << endl;
    }
    LOG << "Vb::Done writing results." << endl;
}

void Vb::CheckCoordMatrixCorrectlyOrdered(const Matrix &coords)
{
    // Only 3D
    assert(coords.Nrows() == 3);

    // Go through each voxel one at a time apart from last
    for (int v = 1; v <= coords.Ncols() - 1; v++)
    {
        // Find difference between current coords and next
        ColumnVector diff = coords.Column(v + 1) - coords.Column(v);

        // Check order
        // +1 = +x, +10 = +y, +100 = +z, -100 = -z+x, etc.
        int d = sign(diff(1)) + 10 * sign(diff(2)) + 100 * sign(diff(3));
        if (d <= 0)
        {
            LOG << "Vb::Found mis-ordered voxels " << v << " and " << v + 1 << ": d=" << d << endl;
            throw FabberInternalError(
                "Coordinate matrix must be in correct order to use adjacency-based priors.");
        }
    }
}

// Binary search for data(index) == num
// Assumes data is sorted ascending!!
// Either returns an index such that data(index) == num
//   or -1 if num is not present in data.
static inline int binarySearch(const ColumnVector &data, int num)
{
    int first = 1, last = data.Nrows();

    while (first <= last)
    {
        int test = (first + last) / 2;

        if (data(test) < num)
        {
            first = test + 1;
        }
        else if (data(test) > num)
        {
            last = test - 1;
        }
        else if (data(test) == num)
        {
            return test;
        }
        else
        {
            assert(false); // logic error!  data wasn't sorted?
        }
    }
    return -1;
}

/**
 * Calculate nearest and second-nearest neighbours for the voxels
 */
void Vb::CalcNeighbours(const Matrix &coords, int spatial_dims, vector<vector<int> > &neighbours, vector<vector<int> > &neighbours2)
{
    int nvoxels = coords.Ncols();
    if (nvoxels == 0)
        return;

    // Voxels must be ordered by increasing z, y and x values respectively
    // otherwise binary search for voxel by offset will not work
    CheckCoordMatrixCorrectlyOrdered(coords);

    // Create a column vector with one entry per voxel.
    ColumnVector offsets(nvoxels);

    // Populate offsets with the offset into the
    // matrix of each voxel. We assume that co-ordinates
    // could be zero but not negative
    int xsize = coords.Row(1).Maximum() + 1;
    int ysize = coords.Row(2).Maximum() + 1;
    for (int v = 1; v <= nvoxels; v++)
    {
        int x = coords(1, v);
        int y = coords(2, v);
        int z = coords(3, v);
        int offset = z * xsize * ysize + y * xsize + x;
        offsets(v) = offset;
    }

    // Delta is a list of offsets to find nearest
    // neighbours in x y and z direction (not diagonally)
    // Of course applying these offsets naively would not
    // always work, e.g. offset of -1 in the x direction
    // will not be a nearest neighbour for the first voxel
    // so need to check for this in subsequent code
    vector<int> delta;
    delta.push_back(1);              // next row
    delta.push_back(-1);             // prev row
    delta.push_back(xsize);          // next column
    delta.push_back(-xsize);         // prev column
    delta.push_back(xsize * ysize);  // next slice
    delta.push_back(-xsize * ysize); // prev slice

    // Don't look for neighbours in all dimensions.
    // For example if spatialDims=2, max_delta=3 so we
    // only look for neighbours in rows and columns
    //
    // However note we still need the full list of 3D deltas for later
    int max_delta = spatial_dims * 2 - 1;

    // Neighbours is a vector of vectors, so each voxel
    // will have an entry which is a vector of its neighbours
    neighbours.resize(nvoxels);

    // Go through each voxel. Note that offsets is indexed from 1 not 0
    // however the offsets themselves (potentially) start at 0.
    for (int vid = 1; vid <= nvoxels; vid++)
    {
        // Get the voxel offset into the matrix
        int pos = int(offsets(vid));

        // Now search for neighbours
        for (int n = 0; n <= max_delta; n++)
        {
            // is there a voxel at this neighbour position?
            // indexed from 1; id == -1 if not found.
            int id = binarySearch(offsets, pos + delta[n]);

            // No such voxel: continue
            if (id < 0)
                continue;

            // Check for wrap-around

            // Don't check for wrap around on final co-ord
            // PREVIOUSLY		if (delta.size() >= n + 2)
            // Changed (fixed)? because if spatialDims != 3 we still need
            // to check for wrap around in y-coordinate FIXME check
            if (n < 4)
            {
                bool ignore = false;
                if (delta[n] > 0)
                {
                    int test = delta[n + 2];
                    if (test > 0)
                        ignore = (pos % test) >= test - delta[n];
                }
                else
                {
                    int test = -delta[n + 2];
                    if (test > 0)
                        ignore = (pos % test) < -delta[n];
                }
                if (ignore)
                {
                    continue;
                }
            }

            // If we get this far, add it to the list
            neighbours.at(vid - 1).push_back(id);
        }
    }

    // Similar algorithm but looking for Neighbours-of-neighbours, excluding self,
    // but including duplicates if there are two routes to get there
    // (diagonally connected)
    neighbours2.resize(nvoxels);

    for (int vid = 1; vid <= nvoxels; vid++)
    {
        // Go through the list of neighbours for each voxel.
        for (unsigned n1 = 0; n1 < neighbours.at(vid - 1).size(); n1++)
        {
            // n1id is the voxel index (not the offset) of the neighbour
            int n1id = neighbours[vid - 1].at(n1);
            int checkNofN = 0;
            // Go through each of it's neighbours. Add each, apart from original voxel
            for (unsigned n2 = 0; n2 < neighbours.at(n1id - 1).size(); n2++)
            {
                int n2id = neighbours[n1id - 1].at(n2);
                if (n2id != vid)
                {
                    neighbours2[vid - 1].push_back(n2id);
                }
                else
                    checkNofN++;
            }

            if (checkNofN != 1)
            {
                throw FabberInternalError("Each of this voxel's neighbours must have "
                                          "this voxel as a neighbour");
            }
        }
    }
}

void SpatialVbThreadContext::Run()
{
    // Pass in some (dummy) data/coords here just in case the model relies upon it
    // use the first voxel values as our dummies FIXME this shouldn't really be
    // necessary, need to find way for model to know about the data beforehand.
    if (nvoxels > 0)
        PassModelData(1);

    vector<Parameter> params;
    m_model->GetParameters(*m_rundata, params);
    vector<Prior *> priors = PriorFactory(*m_rundata).CreatePriors(params);

    // Spatial loop currently uses a global convergence detector FIXME
    // needs to change
    CountingConvergenceDetector conv;
    conv.Initialize(*m_rundata);
    double Fglobal = 1234.5678;
    int maxits = m_rundata->GetIntDefault("max-iterations", 10);

    // MAIN ITERATION LOOP
    do
    {
        LOG << endl << "*** Spatial iteration *** " << (it + 1) << endl;

        // Give an indication of the progress through the voxels;
        m_rundata->Progress(it, maxits);
        double Fprior = 0;

        // ITERATE OVER VOXELS
        for (int v = 1; v <= nvoxels; v++)
        {
            v = v;

            PassModelData(v);

            // The steps below are essentially the same as regular VB, although
            // the code looks different as the per-voxel dists are set up at the
            // start rather than as we go
            try
            {
                Fprior = 0;

                // Apply prior updates for spatial or ARD priors
                for (int k = 0; k < m_num_params; k++)
                {
                    Fprior += priors[k]->ApplyToMVN(&fwd_prior[v - 1], *this);
                }
                if (m_debug)
                    DebugVoxel(v, "Priors set");

                // Ignore voxels where numerical issues have occurred
                if (std::find(ignore_voxels.begin(), ignore_voxels.end(), v)
                    != ignore_voxels.end())
                {
                    LOG << "Ignoring voxel " << v << endl;
                    continue;
                }

                CalculateF(v, "before", Fprior);

                m_noise->UpdateTheta(*noise_post[v - 1], fwd_post[v - 1],
                    fwd_prior[v - 1], m_lin_model[v - 1], m_origdata.Column(v), NULL, 0);
                if (m_debug)
                    DebugVoxel(v, "Theta updated");

                CalculateF(v, "theta", Fprior);
            }
            catch (FabberInternalError &e)
            {
                LOG << "Vb::Internal error for voxel " << v << " at " << m_coords.Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
            catch (NEWMAT::Exception &e)
            {
                LOG << "Vb::NEWMAT exception for voxel " << v << " at " << m_coords.Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
        }

        Fglobal = 0;
        for (int v = 1; v <= nvoxels; v++)
        {
            try {
                // Ignore voxels where numerical issues have occurred
                if (std::find(ignore_voxels.begin(), ignore_voxels.end(), v)
                    != ignore_voxels.end())
                {
                    LOG << "Ignoring voxel " << v << endl;
                    continue;
                }

                PassModelData(v);

                m_noise->UpdateNoise(*noise_post[v - 1], *noise_prior[v - 1],
                    fwd_post[v - 1], m_lin_model[v - 1], m_origdata.Column(v));
                if (m_debug)
                    DebugVoxel(v, "Noise updated");

                CalculateF(v, "noise", Fprior);

                if (!m_locked_linear)
                {
                    m_lin_model[v - 1].ReCentre(fwd_post[v - 1].means);
                    if (m_debug)
                        DebugVoxel(v, "Re-centre");
                }
                
                Fglobal += CalculateF(v, "lin", Fprior);
            }
            catch (FabberInternalError &e)
            {
                LOG << "Vb::Internal error for voxel " << v << " at " << m_coords.Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
            catch (NEWMAT::Exception &e)
            {
                LOG << "Vb::NEWMAT exception for voxel " << v << " at " << m_coords.Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
        }

        ++it;
    } while (!conv.Test(Fglobal));

    // Interesting addition: calculate "coefficient resels" from Penny et al. 2005
    for (int k = 1; k <= m_num_params; k++)
    {
        ColumnVector gamma_vk(nvoxels);
        for (int v = 1; v <= nvoxels; v++)
        {
            gamma_vk(v) = 1
                - fwd_post[v - 1].GetCovariance()(k, k)
                    / fwd_prior[v - 1].GetCovariance()(k, k);
        }
        LOG << "Vb::Coefficient resels per voxel for param " << k << ": "
            << gamma_vk.Sum() / nvoxels << endl;
    }

    for (int v = 1; v <= nvoxels; v++)
    {
        resultMVNs[v - 1]
            = new MVNDist(fwd_post[v - 1], noise_post[v - 1]->OutputAsMVN());
    }
    for (unsigned int i = 0; i < priors.size(); i++)
    {
        delete priors[i];
    }
}

VbThreadContext::VbThreadContext(FabberRunData &rundata, int worker_id, int n_workers, int start_vox, int num_vox)
    : ThreadContext(rundata, worker_id, n_workers, start_vox, num_vox)
{
    // Free energy
    m_saveF = rundata.GetBool("save-free-energy");
    m_printF = rundata.GetBool("print-free-energy");
    if (nvoxels > 0) 
    {
        m_needF = m_conv[0]->UseF() || m_printF || m_saveF;
    }
    else 
    {
        m_needF = m_printF || m_saveF;
    }
}

void VbThreadContext::Run()
{
    vector<Parameter> params;
    m_model->GetParameters(*m_rundata, params);
    vector<Prior *> priors = PriorFactory(*m_rundata).CreatePriors(params);

    // Loop over voxels
    for (int v = 1; v <= nvoxels; v++)
    {
        PassModelData(v);

        v = v;
        it = 0;

        // Save our model parameters in case we need to revert later.
        // Note need to save prior in case ARD is being used
        NoiseParams *const noisePosteriorSave = noise_post[v - 1]->Clone();
        MVNDist fwdPosteriorSave(fwd_post[v - 1]);
        MVNDist fwdPriorSave(fwd_prior[v - 1]);

        // Give an indication of the progress through the voxels;
        m_rundata->Progress(v, nvoxels);
        double F = 1234.5678;

        try
        {
            m_lin_model[v - 1].ReCentre(fwd_post[v - 1].means);
            m_conv[v - 1]->Reset();

            // START the VB updates and run through the relevant iterations (according to the
            // convergence testing)
            do
            {
                double Fprior = 0;

                if (m_conv[v - 1]->NeedRevert()) // revert to previous solution if the convergence
                                                 // detector calls for it
                {
                    *noise_post[v - 1] = *noisePosteriorSave;
                    fwd_post[v - 1] = fwdPosteriorSave;
                    fwd_prior[v - 1] = fwdPriorSave;
                    m_lin_model[v - 1].ReCentre(fwd_post[v - 1].means);
                    if (m_debug)
                        DebugVoxel(v, "Reverted");
                }

                // Save old values if called for
                if (m_conv[v - 1]->NeedSave())
                {
                    *noisePosteriorSave = *noise_post[v - 1]; // copy values, not pointer!
                    fwdPosteriorSave = fwd_post[v - 1];
                    fwdPriorSave = fwd_prior[v - 1];
                }

                for (int k = 0; k < m_num_params; k++)
                {
                    Fprior += priors[k]->ApplyToMVN(&fwd_prior[v - 1], *this);
                }

                if (m_debug)
                    DebugVoxel(v, "Applied priors");

                F = CalculateF(v, "before", Fprior);

                m_noise->UpdateTheta(*noise_post[v - 1], fwd_post[v - 1],
                    fwd_prior[v - 1], m_lin_model[v - 1], m_origdata.Column(v), NULL,
                    m_conv[v - 1]->LMalpha());

                if (m_debug)
                    DebugVoxel(v, "Updated params");

                F = CalculateF(v, "theta", Fprior);

                m_noise->UpdateNoise(*noise_post[v - 1], *noise_prior[v - 1],
                    fwd_post[v - 1], m_lin_model[v - 1], m_origdata.Column(v));

                if (m_debug)
                    DebugVoxel(v, "Updated noise");

                F = CalculateF(v, "phi", Fprior);

                // Linearization update
                // Update the linear model before doing Free energy calculation
                // (and ready for next round of theta and phi updates)
                m_lin_model[v - 1].ReCentre(fwd_post[v - 1].means);

                if (m_debug)
                    DebugVoxel(v, "Re-centered");

                F = CalculateF(v, "lin", Fprior);

                ++it;
            } while (!m_conv[v - 1]->Test(F));

            if (m_debug)
                LOG << "Converged after " << it << " iterations" << endl;

            // Revert to old values at last stage if required
            if (m_conv[v - 1]->NeedRevert())
            {
                *noise_post[v - 1] = *noisePosteriorSave;
                fwd_post[v - 1] = fwdPosteriorSave;
                fwd_prior[v - 1] = fwdPriorSave;
                m_lin_model[v - 1].ReCentre(fwd_post[v - 1].means);
            }

            delete noisePosteriorSave;
        }
        catch (FabberInternalError &e)
        {
            LOG << "Vb::Internal error for voxel " << v << " at " << m_coords.Column(v).t()
                << " : " << e.what() << endl;

            if (m_halt_bad_voxel)
                throw;
        }
        catch (NEWMAT::Exception &e)
        {
            LOG << "Vb::NEWMAT exception for voxel " << v << " at " << m_coords.Column(v).t()
                << " : " << e.what() << endl;

            if (m_halt_bad_voxel)
                throw;
        }

        // now write the results to resultMVNs
        try
        {
            resultMVNs.at(v - 1)
                = new MVNDist(fwd_post[v - 1], noise_post[v - 1]->OutputAsMVN());
                
            if (m_printF)
                resultFs.at(v - 1) = F;
        }
        catch (...)
        {
            // Even that can fail, due to results being singular
            LOG << "Vb::Can't give any sensible answer for this voxel; outputting zero +- "
                   "identity\n";
            MVNDist *tmp = new MVNDist(m_log);
            tmp->SetSize(fwd_post[v - 1].means.Nrows()
                + noise_post[v - 1]->OutputAsMVN().means.Nrows());
            tmp->SetCovariance(IdentityMatrix(tmp->means.Nrows()));
            resultMVNs.at(v - 1) = tmp;
            if (m_printF)
                resultFs.at(v - 1) = F;
        }
    }
    for (unsigned int i = 0; i < priors.size(); i++)
    {
        delete priors[i];
    }
}

/**
 * Calculate free energy. Note that this is currently unused in spatial VB
 */
double VbThreadContext::CalculateF(int v, string label, double Fprior)
{
    double F = 1234.5678;
    if (m_printF)
    {
        F = m_noise->CalcFreeEnergy(*noise_post[v - 1], *noise_prior[v - 1],
            fwd_post[v - 1], fwd_prior[v - 1], m_lin_model[v - 1],
            m_origdata.Column(v));
        F += Fprior;
        resultFs[v - 1] = F;
        if (m_printF)
        {
            LOG << "Vb::F" << label << " = " << F << endl;
        }
    }
    return F;
}

void VbThreadContext::DebugVoxel(int v, const string &where)
{
    LOG << where << " - voxel " << v << " of " << nvoxels << endl;
    LOG << "Prior means: " << endl << fwd_prior[v - 1].means.t();
    LOG << "Prior precisions: " << endl << fwd_prior[v - 1].GetPrecisions();
    LOG << "Posterior means: " << endl << fwd_post[v - 1].means.t();
    LOG << "Noise prior means: " << endl << noise_prior[v - 1]->OutputAsMVN().means.t();
    LOG << "Noise prior precisions: " << endl
        << noise_prior[v - 1]->OutputAsMVN().GetPrecisions();
    LOG << "Centre: " << endl << m_lin_model[v - 1].Centre();
    LOG << "Offset: " << endl << m_lin_model[v - 1].Offset();
    LOG << "Jacobian: " << endl << m_lin_model[v - 1].Jacobian() << endl;
}
