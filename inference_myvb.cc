/*  inference_spatialvb.cc - implementation of VB with spatial priors

 Adrian Groves and Matthew Webster, FMRIB Image Analysis Group

 Copyright (C) 2007-2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_myvb.h"

#include "priors.h"
#include "run_context.h"
#include "convergence.h"
#include "easylog.h"
#include "tools.h"

#include <miscmaths/miscmaths.h>
#include <newmatio.h>

#include <math.h>
#include <iomanip>

using MISCMATHS::sign;

static OptionSpec OPTIONS[] = {
    { "spatial-dims", OPT_INT, "Number of spatial dimensions", OPT_NONREQ, "3" },
    { "spatial-speed", OPT_STR, "Restrict speed of spatial smoothing", OPT_NONREQ,
        "-1" },
    { "distance-measure", OPT_STR, "", OPT_NONREQ, "dist1" },
    { "param-spatial-priors", OPT_STR,
        "Type of spatial priors for each parameter, as a sequence of characters. "
        "N=nonspatial, M=Markov random field, P=Penny, A=ARD",
        OPT_NONREQ, "S+" },
    { "update-spatial-prior-on-first-iteration", OPT_BOOL, "", OPT_NONREQ, "" },
    { "locked-linear-from-mvn", OPT_MVN, "MVN file containing fixed centres for linearization", OPT_NONREQ,
        "" },
    { "" },
};

void Vb::GetOptions(vector<OptionSpec> &opts) const
{
    VariationalBayesInferenceTechnique::GetOptions(opts);

    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

InferenceTechnique *Vb::NewInstance()
{
    return new Vb();
}

string Vb::GetPriorTypesString(FabberRunData &rundata)
{
    string priors_str = rundata.GetStringDefault("param-spatial-priors", "");
    
    // Find out how many prior types are in the string, and what the + character
    // should be interpreted as
    int n_str_params = 0;
    char repeat_type = '-';
    bool no_plus=true;
    for (size_t i=0; i<priors_str.size(); i++) 
    {
        if (priors_str[i] != '+') {
            if (no_plus) repeat_type = priors_str[i];
            n_str_params++;
        }
        else {
            no_plus = false;
        }
    }

    if (n_str_params > m_num_params)
    {
        throw InvalidOptionValue("param-spatial-priors", priors_str, "Too many parameters");
    }

    if (int(priors_str.size()) < m_num_params)
    {
        // Expand '+' char, if present, to give correct number of parameters
        // If there is no +, append with '-', meaning 'model default'
        int deficit = m_num_params - priors_str.size();
        size_t plus_pos = priors_str.find("+");
        if (plus_pos != std::string::npos)
        {
            priors_str.insert(plus_pos, deficit, '+');
        }
        else {
            priors_str.insert(priors_str.end(), deficit, '-');
        }
    }

    // Finally, replace all + chars with identified repeat type    
    std::replace(priors_str.begin(), priors_str.end(), '+', repeat_type);

    assert(int(priors_str.size()) == m_num_params);
    return priors_str;
}
/**
 * How many spatial dimensions have been requested?
 */
static int GetSpatialDims(FabberRunData &args)
{
    int dims = args.GetIntDefault("spatial-dims", 3);

    if (dims < 0 || dims > 3)
    {
        throw InvalidOptionValue("spatial-dims", stringify(dims),
            "Must be 0, 1, 2 or 3");
    }
    else if (dims == 1)
    {
        // WARN_ONCE("spatial-dims=1 is very weird... I hope you're just testing
        // something!");
    }
    else if (dims == 2)
    {
        // WARN_ONCE("spatial-dims=2 doesn't decompose into slices and won't help if
        // you're using the D prior");
    }

    return dims;
}

void Vb::Initialize(FwdModel *fwd_model, FabberRunData &args)
{
    VariationalBayesInferenceTechnique::Initialize(fwd_model, args);

    m_prior_types_str = GetPriorTypesString(args);
    args.Set("param-spatial-priors", m_prior_types_str);
    m_spatial_dims = GetSpatialDims(args);

    // Some unsupported options:

    // Locked linearizations, if requested
    m_locked_linear_file = args.GetStringDefault("locked-linear-from-mvn", "");
    m_locked_linear = (m_locked_linear_file != "");
}

void Vb::SetupPerVoxelDists(FabberRunData &rundata)
{
    // Initialized in voxel loop below (from file or default as required)
    m_noise_post.resize(m_nvoxels, NULL);
    m_noise_prior.resize(m_nvoxels, NULL);
    m_fwd_post.resize(m_nvoxels);

    // Re-centred in voxel loop below
    m_lin_model.resize(m_nvoxels, LinearizedFwdModel(m_model));

    // Static initialization for all voxels currently
    m_fwd_prior.resize(m_nvoxels, MVNDist(m_num_params, m_log));

    // Loaded from file if required, otherwise initialized during calculation
    resultMVNs.resize(m_nvoxels, NULL);

    // Initialized during calculation
    resultFs.resize(m_nvoxels, 9999); // 9999 is a garbage default value

    // Whether to fix the linearization centres (default: false)
    vector<MVNDist *> lockedLinearDists;
    Matrix lockedLinearCentres;
    const int nNoiseParams = initialNoisePrior->OutputAsMVN().GetSize();
    if (m_locked_linear)
    {
        LOG << "Vb::Loading fixed linearization centres from the MVN '"
            << m_locked_linear_file << "'\nNOTE: This does not check if the correct "
                                   "number of parameters is present!\n";
        MVNDist::Load(lockedLinearDists, m_locked_linear_file, rundata, m_log);
        lockedLinearCentres.ReSize(m_num_params, m_nvoxels);
    }

    if (m_continueFromFile != "")
    {
        LOG << "Vb::Continuing from file "
            << m_continueFromFile << endl;
        InitMVNFromFile(m_continueFromFile, rundata, paramFilename);
    }

    for (int v = 1; v <= m_nvoxels; v++)
    {
        if (m_continueFromFile != "")
        {
            m_fwd_post[v - 1] = resultMVNs.at(v - 1)->GetSubmatrix(1, m_num_params);
            assert(m_num_params + nNoiseParams == resultMVNs.at(v - 1)->GetSize());
            m_noise_post[v - 1] = m_noise->NewParams();
            m_noise_post[v - 1]->InputFromMVN(resultMVNs.at(v - 1)->GetSubmatrix(
                m_num_params + 1, m_num_params + nNoiseParams));
        }
        else
        {
            // Set the initial posterior for model params. Model
            // may want the voxel data in order to do this
            PassModelData(v); 
            m_model->GetInitialPosterior(m_fwd_post[v - 1]);
            // Set initial noise posterior
            m_noise_post[v - 1] = initialNoisePosterior->Clone();
        }

        if (m_locked_linear)
        {
            lockedLinearCentres.Column(v) = lockedLinearDists.at(v - 1)->means.Rows(1, m_num_params);
            m_lin_model[v - 1].ReCentre(lockedLinearCentres.Column(v));
        }
        else
        {
            //LOG  << "Initial re-centering: " << m_fwd_post[v - 1].means.t() << endl;
            m_lin_model[v - 1].ReCentre(m_fwd_post[v - 1].means);
        }

        m_noise_prior[v - 1] = initialNoisePrior->Clone();
        m_noise->Precalculate(*m_noise_post[v - 1], *m_noise_prior[v - 1],
            m_origdata->Column(v));
    }
}

void Vb::IgnoreVoxel(int v) 
{
    LOG << "Vb::IgnoreVoxel This voxel will be ignored in further updates" << endl;
    
    m_ignore_voxels.push_back(v);
    
    // Remove voxel from lists of neighbours of other voxels.
    // We identify affected voxels by looking in the neighbour
    // lists for the bad voxel, because any voxel which has 
    // the bad voxel as a neighbour will be a neighbour of the
    // bad voxel
    vector<int> nn = m_neighbours[v-1];
    for (vector<int>::iterator i=nn.begin(); i!=nn.end(); ++i) {
        // Reference to list of neighbours of some other voxel which
        // has the bad voxel as a neighbour
        vector<int> &n2 = m_neighbours[*i-1];

        n2.erase(std::remove(n2.begin(), n2.end(), v), n2.end());
    }

    // Same for next-nearest-neighbours
    nn = m_neighbours2[v-1];
    for (vector<int>::iterator i=nn.begin(); i!=nn.end(); ++i) {
        // Reference to list of neighbours of some other voxel which
        // has the bad voxel as a neighbour
        vector<int> &n2 = m_neighbours2[*i-1];

        n2.erase(std::remove(n2.begin(), n2.end(), v), n2.end());
    }
}

/**
 * Calculate free energy. Note that this is currently unused in spatial VB
 */
double Vb::CalculateF(int v, string label, double Fprior)
{
    double F = 1234.5678;
    if (m_needF)
    {
        F = m_noise->CalcFreeEnergy(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                    m_fwd_post[v - 1], m_fwd_prior[v - 1],
                    m_lin_model[v - 1], m_origdata->Column(v));
        F += Fprior;
        resultFs[v-1] = F;
        if (m_printF)
        {
            LOG << "Vb::F" << label << " = " << F << endl;
        }
    }
    return F;
}

// M = Markov random field - normally used
// P = Alternative to M (Penny prior?)
// N = non-spatial prior (model default)
// A = ARD prior
// I = image prior
//
// m/p are variations on M/P with different edge behaviour (Dirichlet BCs)
void Vb::DoCalculations(FabberRunData &rundata)
{    
    // extract data (and the coords) from rundata for the (first) VB run
    // Rows are volumes
    // Columns are (time) series
    // num Rows is size of (time) series
    // num Cols is size of volumes
    m_origdata = &rundata.GetMainVoxelData();
    m_coords = &rundata.GetVoxelCoords();
    m_suppdata = &rundata.GetVoxelSuppData();
    m_nvoxels = m_origdata->Ncols();

    // pass in some (dummy) data/coords here just in case the model relies upon it
    // use the first voxel values as our dummies FIXME this shouldn't really be
    // necessary, need to find way for model to know about the data beforehand.
    if (m_nvoxels > 0)
        PassModelData(1);

    // Only call DoCalculations once
    assert(resultMVNs.empty());
    assert(resultFs.empty());

    // FIXME hack
    vector<Parameter> params;
    m_model->GetParameters(rundata, params);
    SetupPerVoxelDists(rundata);

    if (rundata.GetBool("output-only"))
    {
        // Do no calculations - now we have set resultMVNs we can finish
        LOG << "Vb::DoCalculations output-only set - not performing any calculations" << endl;
    }
    else if (rundata.GetString("method") == "vb") {
        DoCalculationsVoxelwise(rundata);
    }
    else {
        DoCalculationsSpatial(rundata);        
    }

    if (!m_needF)
    {
        // clearing resultFs here should prevent an F image from being saved.
        resultFs.clear();
    }

    // Delete stuff (avoid memory leaks)
    for (int v = 1; v <= m_nvoxels; v++)
    {
        delete m_noise_post[v - 1];
        delete m_noise_prior[v - 1];
    }
}

void Vb::DoCalculationsVoxelwise(FabberRunData &rundata)
{
    vector<Parameter> params;
    m_model->GetParameters(rundata, params);
    vector<Prior *> priors = PriorFactory(rundata).CreatePriors(params);

    RunContext ctx(m_nvoxels, m_fwd_post, m_neighbours, m_neighbours2);
    m_conv->Reset();
    
    // Loop over voxels
    for (int v = 1; v <= m_nvoxels; v++)
    {
        PassModelData(v);

        ctx.v = v;
        ctx.it = 0; 

        // Save our model parameters in case we need to revert later.
        // Note need to save prior in case ARD is being used
        NoiseParams *const noisePosteriorSave = m_noise_post[v - 1]->Clone();
        MVNDist fwdPosteriorSave(m_fwd_post[v - 1]);
        MVNDist fwdPriorSave(m_fwd_prior[v - 1]);

        // Give an indication of the progress through the voxels;
        rundata.Progress(v, m_nvoxels);
        double F = 1234.5678;

        try
        {
            m_lin_model[v-1].ReCentre(m_fwd_post[v - 1].means);
            m_noise->Precalculate(*m_noise_post[v - 1], *m_noise_prior[v - 1], m_origdata->Column(v));
            m_conv->Reset();

            // START the VB updates and run through the relevant iterations (according to the convergence testing)
            do
            {
                double Fprior = 0;
                
                if (m_conv->NeedRevert()) //revert to previous solution if the convergence detector calls for it
                {
                    *m_noise_post[v - 1] = *noisePosteriorSave;
                    m_fwd_post[v - 1] = fwdPosteriorSave;
                    m_fwd_prior[v-1] = fwdPriorSave;
                    m_lin_model[v-1].ReCentre(m_fwd_post[v - 1].means);
                }

                // Save old values if called for
                if (m_conv->NeedSave())
                {
                    *noisePosteriorSave = *m_noise_post[v - 1]; // copy values, not pointer!
                    fwdPosteriorSave = m_fwd_post[v - 1];
                    fwdPriorSave = m_fwd_prior[v-1];
                }
                
                for (int k=0; k<m_num_params; k++) {
                    Fprior += priors[k]->ApplyToMVN(&m_fwd_prior[v-1], ctx);
                }

                F = CalculateF(v, "before", Fprior);

                m_noise->UpdateTheta(*m_noise_post[v - 1], m_fwd_post[v - 1], m_fwd_prior[v-1], m_lin_model[v-1], m_origdata->Column(v), NULL,
                    m_conv->LMalpha());

                F = CalculateF(v, "theta", Fprior);

                m_noise->UpdateNoise(*m_noise_post[v - 1], *m_noise_prior[v - 1], m_fwd_post[v - 1], m_lin_model[v-1], m_origdata->Column(v));
                
                F = CalculateF(v, "phi", Fprior);
                
                // Linearization update
                // Update the linear model before doing Free energy calculation 
                // (and ready for next round of theta and phi updates)
                m_lin_model[v-1].ReCentre(m_fwd_post[v - 1].means);

                F = CalculateF(v, "lin", Fprior);
                
                ++ctx.it;
            } while (!m_conv->Test(F));

            // Revert to old values at last stage if required
            if (m_conv->NeedRevert())
            {
                *m_noise_post[v - 1] = *noisePosteriorSave;
                m_fwd_post[v - 1] = fwdPosteriorSave;
                m_fwd_prior[v-1] = fwdPriorSave;
                m_lin_model[v-1].ReCentre(m_fwd_post[v - 1].means);
            }
        }
        catch (FabberInternalError &e) {
            LOG << "Vb::Internal error for voxel " <<  v 
                << " at " << m_coords->Column(v).t() << " : " << e.what() << endl;
            
            if (m_halt_bad_voxel) throw;
            else IgnoreVoxel(v);
        }
        catch (NEWMAT::Exception &e) {

            LOG << "Vb::NEWMAT exception for voxel " <<  v 
                << " at " << m_coords->Column(v).t() << " : " << e.what() << endl;
            
            if (m_halt_bad_voxel) throw;
            else IgnoreVoxel(v);
        }

        // now write the results to resultMVNs
        try
        {
            resultMVNs.at(v - 1) = new MVNDist(m_fwd_post[v-1], m_noise_post[v-1]->OutputAsMVN());
            if (m_needF)
            {
                resultFs.at(v - 1) = F;
            }
        }
        catch (...)
        {
            // Even that can fail, due to results being singular
            LOG
                << "Vb::Can't give any sensible answer for this voxel; outputting zero +- identity\n";
            MVNDist *tmp = new MVNDist(m_log);
            tmp->SetSize(m_fwd_post[v-1].means.Nrows() + m_noise_post[v-1]->OutputAsMVN().means.Nrows());
            tmp->SetCovariance(IdentityMatrix(tmp->means.Nrows()));
            resultMVNs.at(v - 1) = tmp;

            if (m_needF)
                resultFs.at(v - 1) = F;
        }
    }
}

void Vb::DoCalculationsSpatial(FabberRunData &rundata)
{
   int maxits = convertTo<int>(rundata.GetStringDefault("max-iterations", "10"));
   if (m_debug) LOG << setprecision(17);

    // pass in some (dummy) data/coords here just in case the model relies upon it
    // use the first voxel values as our dummies FIXME this shouldn't really be
    // necessary, need to find way for model to know about the data beforehand.
    if (m_nvoxels > 0)
        PassModelData(1);

    // Make the neighbours[] lists if required
    //if (m_prior_types_str.find_first_of("mMpP") != string::npos)
    if (true) // FIXME
    {
        CalcNeighbours(*m_coords);
    }

    vector<Parameter> params;
    m_model->GetParameters(rundata, params);
    vector<Prior *> priors = PriorFactory(rundata).CreatePriors(params);
    
    RunContext ctx(m_nvoxels, m_fwd_post, m_neighbours, m_neighbours2);
    m_conv->Reset();
    
    double Fglobal=1234.5678;

    // MAIN ITERATION LOOP
    do
    {
        LOG << endl << "*** Spatial iteration *** " << (ctx.it+1) << endl;

        // Give an indication of the progress through the voxels;
        rundata.Progress(ctx.it, maxits);
        double Fprior=0;

        // ITERATE OVER VOXELS
        for (int v = 1; v <= m_nvoxels; v++)
        {
            ctx.v = v;

            PassModelData(v);

            // The steps below are essentially the same as regular VB, although
            // the code looks different as the per-voxel dists are set up at the
            // start rather than as we go
            try {
                Fprior=0;

                // Apply prior updates for spatial or ARD priors
                for (int k=0; k<m_num_params; k++) {
                    Fprior += priors[k]->ApplyToMVN(&m_fwd_prior[v-1], ctx);
                }
                if (m_debug) {
                    LOG << "Voxel " << v << " of " << m_nvoxels << endl;
                    LOG << "Prior means: " << m_fwd_prior[v-1].means.t();
                    LOG << "Prior precisions: " << m_fwd_prior[v-1].GetPrecisions();
                }

                // Ignore voxels where numerical issues have occurred
                if (std::find(m_ignore_voxels.begin(), m_ignore_voxels.end(), v) != m_ignore_voxels.end()) continue;

                CalculateF(v, "before", Fprior);

                m_noise->UpdateTheta(*m_noise_post[v - 1], m_fwd_post[v - 1],
                    m_fwd_prior[v - 1], m_lin_model[v - 1],
                    m_origdata->Column(v), NULL, 0);
                if (m_debug) {
                    LOG << "Post means: " << m_fwd_post[v-1].means.t();
                    LOG << "Post precisions: " << m_fwd_post[v-1].GetPrecisions();
                }
                CalculateF(v, "theta", Fprior);
            }
            catch (FabberInternalError &e) {
                LOG << "Vb::Internal error for voxel " <<  v 
                    << " at " << m_coords->Column(v).t() << " : " << e.what() << endl;
                
                if (m_halt_bad_voxel) throw;
                else IgnoreVoxel(v);
            }
            catch (NEWMAT::Exception &e) {

                LOG << "Vb::NEWMAT exception for voxel " <<  v 
                    << " at " << m_coords->Column(v).t() << " : " << e.what() << endl;
                
                if (m_halt_bad_voxel) throw;
                else IgnoreVoxel(v);
            }
        }

        Fglobal=0;
        for (int v = 1; v <= m_nvoxels; v++) 
        {
            PassModelData(v);

            m_noise->UpdateNoise(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                m_fwd_post[v - 1], m_lin_model[v - 1],
                m_origdata->Column(v));

            CalculateF(v, "noise", Fprior);

            if (!m_locked_linear)
                m_lin_model[v - 1].ReCentre(m_fwd_post[v - 1].means);

            Fglobal += CalculateF(v, "lin", Fprior);
        }

        ++ctx.it;
    } while (!m_conv->Test(Fglobal));

    // Interesting addition: calculate "coefficient resels" from Penny et al. 2005
    for (int k = 1; k <= m_num_params; k++)
    {
        ColumnVector gamma_vk(m_nvoxels);
        for (int v = 1; v <= m_nvoxels; v++)
        {
            gamma_vk(v) = 1 - m_fwd_post[v - 1].GetCovariance()(k, k) / m_fwd_prior[v - 1].GetCovariance()(k, k);
        }
        LOG << "Vb::Coefficient resels per voxel for param "
            << k << ": " << gamma_vk.Sum() / m_nvoxels << endl;
    }

    for (int v = 1; v <= m_nvoxels; v++)
    {
        resultMVNs[v - 1] = new MVNDist(m_fwd_post[v - 1], m_noise_post[v - 1]->OutputAsMVN());
    }
}

void Vb::CheckCoordMatrixCorrectlyOrdered(const Matrix &coords)
{
    // Only 3D
    assert(coords.Nrows() == 3);

    // Voxels are stored one per column, each column is the x/y/z coords
    const int m_nvoxels = coords.Ncols();

    // Go through each voxel one at a time apart from last
    for (int v = 1; v <= m_nvoxels - 1; v++)
    {
        // Find difference between current coords and next
        ColumnVector diff = coords.Column(v + 1) - coords.Column(v);

        // Check order
        // +1 = +x, +10 = +y, +100 = +z, -100 = -z+x, etc.
        int d = sign(diff(1)) + 10 * sign(diff(2)) + 100 * sign(diff(3));
        if (d <= 0)
        {
            LOG << "Vb::Found mis-ordered voxels " << v << " and " << v + 1 << ": d=" << d << endl;
            throw FabberInternalError("Coordinate matrix must be in correct order to use adjacency-based priors.");
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
void Vb::CalcNeighbours(const Matrix &coords)
{
    const int nVoxels = coords.Ncols();
    if (nVoxels == 0)
        return;

    // Voxels must be ordered by increasing z, y and x values respectively
    // otherwise binary search for voxel by offset will not work
    CheckCoordMatrixCorrectlyOrdered(coords);

    // Create a column vector with one entry per voxel.
    ColumnVector offsets(nVoxels);

    // Populate offsets with the offset into the
    // matrix of each voxel. We assume that co-ordinates
    // could be zero but not negative
    int xsize = coords.Row(1).Maximum() + 1;
    int ysize = coords.Row(2).Maximum() + 1;
    for (int v = 1; v <= nVoxels; v++)
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
    int max_delta = m_spatial_dims * 2 - 1;

    // Neighbours is a vector of vectors, so each voxel
    // will have an entry which is a vector of its neighbours
    m_neighbours.resize(nVoxels);

    // Go through each voxel. Note that offsets is indexed from 1 not 0
    // however the offsets themselves (potentially) start at 0.
    for (int vid = 1; vid <= nVoxels; vid++)
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
            m_neighbours.at(vid - 1).push_back(id);
        }
    }

    // Similar algorithm but looking for Neighbours-of-neighbours, excluding self,
    // but including duplicates if there are two routes to get there
    // (diagonally connected) 
    m_neighbours2.resize(nVoxels);

    for (int vid = 1; vid <= nVoxels; vid++)
    {
        // Go through the list of neighbours for each voxel.
        for (unsigned n1 = 0; n1 < m_neighbours.at(vid - 1).size(); n1++)
        {
            // n1id is the voxel index (not the offset) of the neighbour
            int n1id = m_neighbours[vid - 1].at(n1);
            int checkNofN = 0;
            // Go through each of it's neighbours. Add each, apart from original voxel
            for (unsigned n2 = 0; n2 < m_neighbours.at(n1id - 1).size(); n2++)
            {
                int n2id = m_neighbours[n1id - 1].at(n2);
                if (n2id != vid)
                {
                    m_neighbours2[vid - 1].push_back(n2id);
                }
                else
                    checkNofN++;
            }

            if (checkNofN != 1)
            {
                throw FabberInternalError(
                    "Each of this voxel's neighbours must have "
                    "this voxel as a neighbour");
            }
        }
    }
}

void Vb::SaveResults(FabberRunData &rundata) const
{
    LOG << "Vb::Preparing to save results..." << endl;
    int nVoxels = resultMVNs.size();

    // Save the resultMVNs
    if (rundata.GetBool("save-mvn"))
    {
        MVNDist::Save(resultMVNs, "finalMVN", rundata);
    }

    // Create individual files for each parameter's mean and Z-stat
    vector<string> paramNames;
    m_model->NameParams(paramNames);

    if (rundata.GetBool("save-mean") | rundata.GetBool("save-std") | rundata.GetBool("save-zstat"))
    {
        LOG << "Vb::Writing means..." << endl;
        for (unsigned i = 1; i <= paramNames.size(); i++)
        {
            Matrix paramMean, paramZstat, paramStd;
            paramMean.ReSize(1, nVoxels);
            paramZstat.ReSize(1, nVoxels);
            paramStd.ReSize(1, nVoxels);
            
            for (int vox = 1; vox <= nVoxels; vox++)
            {
                MVNDist result = *resultMVNs[vox - 1];
                m_model->ToModel(result);
                paramMean(1, vox) = result.means(i);
                double std = sqrt(result.GetCovariance()(i, i));
                paramZstat(1, vox) = paramMean(1, vox) / std;
                paramStd(1, vox) = std;
            }

            if (rundata.GetBool("save-mean"))
                rundata.SaveVoxelData("mean_" + paramNames.at(i - 1), paramMean);
            if (rundata.GetBool("save-zstat"))
                rundata.SaveVoxelData("zstat_" + paramNames.at(i - 1), paramZstat);
            if (rundata.GetBool("save-std"))
                rundata.SaveVoxelData("std_" + paramNames.at(i - 1), paramStd);
        }
    }

    // Produce the model fit and residual volume series
    bool saveModelFit = rundata.GetBool("save-model-fit");
    bool saveResiduals = rundata.GetBool("save-residuals");
    if (saveModelFit || saveResiduals)
    {
        LOG << "Vb::Writing model fit/residuals..." << endl;

        Matrix modelFit, residuals, datamtx, coords;
        datamtx = rundata.GetMainVoxelData(); // it is just possible that the model needs the data in its calculations
        coords = rundata.GetVoxelCoords();
        modelFit.ReSize(datamtx.Nrows(), nVoxels);
        ColumnVector tmp;
        for (int vox = 1; vox <= nVoxels; vox++)
        {
            // pass in stuff that the model might need
            ColumnVector y = datamtx.Column(vox);
            ColumnVector vcoords = coords.Column(vox);
            m_model->PassData(y, vcoords);

            // do the evaluation
            m_model->EvaluateFabber(resultMVNs.at(vox - 1)->means.Rows(1, m_num_params), tmp);
            modelFit.Column(vox) = tmp;
        }

        if (saveResiduals)
        {
            residuals = datamtx - modelFit;
            rundata.SaveVoxelData("residuals", residuals);
        }
        if (saveModelFit)
        {
            rundata.SaveVoxelData("modelfit", modelFit);
        }
    }

    LOG << "Vb::Done writing results." << endl;
}
