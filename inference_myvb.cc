/*  inference_spatialvb.cc - implementation of VB with spatial priors

 Adrian Groves and Matthew Webster, FMRIB Image Analysis Group

 Copyright (C) 2007-2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_myvb.h"

#include "priors.h"
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

static string GetPriorTypesStr(vector<PriorType> &priors)
{
    string types;
    for (size_t i = 0; i < priors.size(); i++)
    {
        types += priors[i].m_type;
    }
    return types + "\0";
}

void Vb::Initialize(FwdModel *fwd_model,
    FabberRunData &args)
{
    // In spatial VB we want to default to spatial priors so we set this parameter
    // if it's not already set by the user. Note that we are doing this before
    // initializing the VB method so it acts as a default
    args.Set("default-prior-type", "N");
    VariationalBayesInferenceTechnique::Initialize(fwd_model, args);

    m_prior_types_str = GetPriorTypesStr(m_prior_types);
    m_spatial_dims = GetSpatialDims(args);

    // Some unsupported options:

    // Locked linearizations, if requested
    m_locked_linear_file = args.GetStringDefault("locked-linear-from-mvn", "");
    m_locked_linear = (m_locked_linear_file != "");
}

void Vb::SetupPerVoxelDists(FabberRunData &allData)
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
        MVNDist::Load(lockedLinearDists, m_locked_linear_file, allData, m_log);
        lockedLinearCentres.ReSize(m_num_params, m_nvoxels);
    }

    if (m_continueFromFile != "")
    {
        LOG << "SpatialVbInferenceTechnique::Continuing from file "
            << m_continueFromFile << endl;
        InitMVNFromFile(m_continueFromFile, allData, paramFilename);
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
            m_fwd_post[v - 1] = *initialFwdPosterior;
            PassModelData(v); // May be required for model to initialize posterior
            m_model->InitParams(m_fwd_post[v - 1]);
            m_noise_post[v - 1] = initialNoisePosterior->Clone();
        }

        if (m_locked_linear)
        {
            lockedLinearCentres.Column(v) = lockedLinearDists.at(v - 1)->means.Rows(1, m_num_params);
            m_lin_model[v - 1].ReCentre(lockedLinearCentres.Column(v));
        }
        else
        {
            m_lin_model[v - 1].ReCentre(m_fwd_post[v - 1].means);
        }

        m_noise_prior[v - 1] = initialNoisePrior->Clone();
        m_noise->Precalculate(*m_noise_post[v - 1], *m_noise_prior[v - 1],
            m_origdata->Column(v));
    }
}

void Vb::IgnoreVoxel(int v) 
{
    LOG << "This voxel will be ignored in further updates" << endl;
    
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
void Vb::CalculateF(int v, string label, double Fprior)
{
    if (m_needF)
    {
        double F = m_noise->CalcFreeEnergy(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                    m_fwd_post[v - 1], m_fwd_prior[v - 1],
                    m_lin_model[v - 1], m_origdata->Column(v));
        F += Fprior;
        resultFs[v-1] = F;
        if (m_printF)
        {
            LOG << "SpatialVbInferenceTechnique::F" << label << " = " << F << endl;
        }
    }
}

// M = Markov random field - normally used
// P = Alternative to M (Penny prior?)
// N = non-spatial prior (model default)
// A = ARD prior
// I = image prior
//
// m/p are variations on M/P with different edge behaviour (Dirichlet BCs)
// 
void Vb::DoCalculations(FabberRunData &allData)
{
    // extract data (and the coords) from allData for the (first) VB run
    // Rows are volumes
    // Columns are (time) series
    // num Rows is size of (time) series
    // num Cols is size of volumes
    m_origdata = &allData.GetMainVoxelData();
    m_coords = &allData.GetVoxelCoords();
    m_suppdata = &allData.GetVoxelSuppData();
    m_nvoxels = m_origdata->Ncols();
    int maxits = convertTo<int>(allData.GetStringDefault("max-iterations", "10"));

    // pass in some (dummy) data/coords here just in case the model relies upon it
    // use the first voxel values as our dummies FIXME this shouldn't really be
    // necessary, need to find way for model to know about the data beforehand.
    if (m_nvoxels > 0)
        PassModelData(1);

    // Only call DoCalculations once
    assert(resultMVNs.empty());
    assert(resultFs.empty());

    // Make the neighbours[] lists if required
    if (m_prior_types_str.find_first_of("mMpP") != string::npos)
    {
        CalcNeighbours(*m_coords);
    }

    SetupPerVoxelDists(allData);
    if (allData.GetBool("output-only"))
    {
        // Do no calculations - now we have set resultMVNs we can finish
        return;
    }

    PriorContext ctx(m_nvoxels, m_fwd_post, m_neighbours, m_neighbours2);
    PriorFactory pfac(*m_model, allData);
    vector<Prior *> priors = pfac.CreatePriors();

    // FIXME can't calculate free energy with spatial VB yet 
    // This value never changes
    const double globalF = 1234.5678;

    m_conv->Reset();
    
    // MAIN ITERATION LOOP
    do
    {
        // Give an indication of the progress through the voxels;
        allData.Progress(ctx.it, maxits);
        double Fprior=0;

        // ITERATE OVER VOXELS
        for (int v = 1; v <= m_nvoxels; v++)
        {
            ctx.v = v;

            // Ignore voxels where numerical issues have occurred
            if (std::find(m_ignore_voxels.begin(), m_ignore_voxels.end(), v) != m_ignore_voxels.end()) continue;

            PassModelData(v);

            Fprior=0;
            for (int k=0; k<m_num_params; k++) {
                Fprior += priors[k]->ApplyToMVN(&m_fwd_prior[v-1], ctx);
            }

            // The steps below are essentially the same as regular VB, although
            // the code looks different as the per-voxel dists are set up at the
            // start rather than as we go

            try {
                CalculateF(v, "before", Fprior);

                m_noise->UpdateTheta(*m_noise_post[v - 1], m_fwd_post[v - 1],
                    m_fwd_prior[v - 1], m_lin_model[v - 1],
                    m_origdata->Column(v), NULL, 0);

                CalculateF(v, "theta", Fprior);
            }
            catch (FabberInternalError &e) {
                LOG << "SpatialVbInferenceTechnique::Internal error for voxel " <<  v 
                    << " at " << m_coords->Column(v).t() << " : " << e.what() << endl;
                
                if (m_halt_bad_voxel) throw;
                else IgnoreVoxel(v);
            }
            catch (NEWMAT::Exception &e) {

                LOG << "SpatialVbInferenceTechnique::NEWMAT exception for voxel " <<  v 
                    << " at " << m_coords->Column(v).t() << " : " << e.what() << endl;
                
                if (m_halt_bad_voxel) throw;
                else IgnoreVoxel(v);
            }
        }

        for (int v = 1; v <= m_nvoxels; v++) 
        {
            PassModelData(v);

            m_noise->UpdateNoise(*m_noise_post[v - 1], *m_noise_prior[v - 1],
                m_fwd_post[v - 1], m_lin_model[v - 1],
                m_origdata->Column(v));

            CalculateF(v, "noise", Fprior);

            // MOVED HERE on Michael's advice -- 2007-11-23
            if (!m_locked_linear)
                m_lin_model[v - 1].ReCentre(m_fwd_post[v - 1].means);

            CalculateF(v, "lin", Fprior);
        }

        ++ctx.it;
    } while (!m_conv->Test(globalF));

    // Interesting addition: calculate "coefficient resels" from Penny et al. 2005
    for (int k = 1; k <= m_num_params; k++)
    {
        ColumnVector gamma_vk(m_nvoxels); // might be handy
        ColumnVector gamma_vk_eo(
            m_nvoxels); // slightly different calculation (differs if using EO)
        gamma_vk_eo = -999;
        for (int v = 1; v <= m_nvoxels; v++)
        {
            gamma_vk(v) = 1 - m_fwd_post[v - 1].GetCovariance()(k, k) / m_fwd_prior[v - 1].GetCovariance()(k, k);
        }
        LOG << "Vb::Coefficient resels per voxel for param "
            << k << ": " << gamma_vk.Sum() / m_nvoxels << " (vb) or "
            << gamma_vk_eo.Sum() / m_nvoxels << " (eo)" << endl;
    }

    for (int v = 1; v <= m_nvoxels; v++)
    {
        resultMVNs[v - 1] = new MVNDist(m_fwd_post[v - 1], m_noise_post[v - 1]->OutputAsMVN());
    }

    // resultFs are stored as we go along.
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
            LOG << "Found mis-ordered voxels " << v << " and " << v + 1 << ": d=" << d << endl;
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
