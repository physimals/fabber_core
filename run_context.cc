/**
 * run_context.cc
 *
 * ThreadContext is currently a pure structure. May wish to add methods later
 *
 * Copyright (C) 2007-2017 University of Oxford
 */

#include "run_context.h"

#include <miscmaths/miscmaths.h>
#include <newmat.h>

#include <fstream>

using std::vector;
using std::string;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;
using NEWMAT::ColumnVector;
using MISCMATHS::sign;

void ThreadContext::InitializeNoiseFromParam(FabberRunData &rundata, NoiseParams *dist, string param_key)
{
    string filename = rundata.GetStringDefault(param_key, "modeldefault");
    if (filename != "modeldefault")
    {
        // FIXME should there be checking of size here
        LOG << "VbInferenceTechnique::Loading " << param_key << " distribution from " << filename
            << endl;
        dist->InputFromMVN(MVNDist(filename, m_log));
    }
}

void ThreadContext::InitMVNFromFile(FabberRunData &rundata)
{
    LOG << "InferenceTechnique::Loading previous output from MVN" << endl;
    Matrix mvn_data = GetVoxelData(rundata, "continue-from-mvn");
    MVNDist::Load(resultMVNs, mvn_data, m_log);
}

void ThreadContext::PassModelData(int v)
{
    // Pass in data, coords and supplemental data for this voxel
    ColumnVector data = m_origdata.Column(v);
    ColumnVector vcoords = m_coords.Column(v);
    if (m_suppdata.Ncols() > 0)
    {
        m_model->PassData(v, data, vcoords, m_suppdata.Column(v));
    }
    else
    {
        m_model->PassData(v, data, vcoords);
    }
}

Matrix ThreadContext::GetVoxelData(FabberRunData &rundata, std::string name)
{
    try 
    {
        Matrix data = rundata.GetVoxelData(name);
        return data.Columns(start_voxel, start_voxel+nvoxels-1);
    }
    catch (DataNotFound &e) 
    {
        Matrix empty;
        return empty;
    }
}

ThreadContext::ThreadContext(FabberRunData &rundata, int worker_id, int num_workers, int start_vox, int num_vox)
    : m_rundata(&rundata)
    , m_id(worker_id)
    , m_num_workers(num_workers)
    , it(0)
    , v(1)
    , start_voxel(1)
    , nvoxels(0)
    , m_num_params(0)
    , m_noise_params(0)
{
    m_log = rundata.GetLogger();
    m_debug = rundata.GetBool("debug");
    m_halt_bad_voxel = !rundata.GetBool("allow-bad-voxels");
    
    if (num_vox < 0) 
    {
        num_vox = rundata.GetMainVoxelData().Ncols();
    }
    start_voxel = start_vox;
    nvoxels = num_vox;

    //cerr << "ThreadContext: " << start_vox << ", " << num_vox << endl;
    m_origdata = rundata.GetMainVoxelData().Columns(start_vox, start_vox+num_vox-1);
    m_coords = rundata.GetVoxelCoords().Columns(start_vox, start_vox+num_vox-1);
    m_suppdata = GetVoxelData(rundata, "suppdata");

    // Read masked time points option and log if any have been specified
    m_masked_tpoints = rundata.GetIntList("mt", 1);
    if (m_masked_tpoints.size() > 0)
    {
        LOG << "ThreadContext::Masking " << m_masked_tpoints.size() << " time points: ";
        for (unsigned int i = 0; i < m_masked_tpoints.size(); i++)
        {
            LOG << m_masked_tpoints[i] << " ";
        }
        LOG << endl;
    }

    m_model = FwdModel::NewFromName(rundata.GetString("model"));
    // For backwards compatibility - model may not have called superclass initialize
    m_model->SetLogger(m_log);
    m_model->Initialize(rundata);
    vector<Parameter> params;
    m_model->GetParameters(rundata, params);
    m_num_params = params.size();

    // Get noise model.
    m_noise = NoiseModel::NewFromName(rundata.GetStringDefault("noise", "white"));
    m_noise->Initialize(rundata);
    if (rundata.GetStringDefault("method", "vb") != "nlls") 
    {
        m_noise_params = m_noise->NumParams();
        //cerr << "noise has " << m_noise_params << " params" << endl;
    }

    // Initialezed in voxel loop below (from file or default as required)
    noise_post.resize(nvoxels, NULL);
    noise_prior.resize(nvoxels, NULL);
    fwd_post.resize(nvoxels);

    // Re-centred in voxel loop below
    m_lin_model.resize(nvoxels, LinearizedFwdModel(m_model));

    // Initialized in voxel loop below
    m_conv.resize(nvoxels, NULL);
    string conv_name = rundata.GetStringDefault("convergence", "maxits");

    // Model prior is updated during main voxel loop
    fwd_prior.resize(nvoxels, MVNDist(m_num_params, m_log));

    // Loaded from file if required, otherwise initialized during calculation
    resultMVNs.resize(nvoxels, NULL);

    // Initialized during calculation
    resultFs.resize(nvoxels, 9999); // 9999 is a garbage default value

    // Whether to fix the linearization centres (default: false)
    string locked_linear = rundata.GetStringDefault("locked-linear-from-mvn", "");
    m_locked_linear = locked_linear != "";
    vector<MVNDist *> lockedLinearDists;
    Matrix lockedLinearCentres;
    if (m_locked_linear)
    {
        string file = rundata.GetString("locked-linear-from-mvn");
        LOG << "Vb::Loading fixed linearization centres from the MVN '" << file
            << "'\nNOTE: This does not check if the correct "
               "number of parameters is present!\n";
        MVNDist::Load(lockedLinearDists, file, rundata, m_log);
        lockedLinearCentres.ReSize(m_num_params, nvoxels);
    }

    // If we are resuming from a previous run, there will be data containing a per-voxel
    // distribution of the model parameters, and noise as well.
    bool continueFromMvn = false;
    try {
        rundata.GetVoxelData("continue-from-mvn");
        continueFromMvn = true;
    }
    catch(DataNotFound &e) {
        // no worries
    }

    if (continueFromMvn)
    {
        LOG << "Vb::Continuing from MVN" << endl;
        InitMVNFromFile(rundata);
    }

    // Initial noise distributions
    auto_ptr<NoiseParams> initialNoisePrior(m_noise->NewParams());
    auto_ptr<NoiseParams> initialNoisePosterior(m_noise->NewParams());
    m_noise->HardcodedInitialDists(*initialNoisePrior, *initialNoisePosterior);
    InitializeNoiseFromParam(rundata, initialNoisePrior.get(), "noise-initial-prior");
    InitializeNoiseFromParam(rundata, initialNoisePosterior.get(), "noise-initial-posterior");

    for (int v = 1; v <= nvoxels; v++)
    {
        if (continueFromMvn)
        {
            fwd_post[v - 1] = resultMVNs.at(v - 1)->GetSubmatrix(1, m_num_params);
            assert(m_num_params + m_noise_params == resultMVNs.at(v - 1)->GetSize());
            noise_post[v - 1] = m_noise->NewParams();
            noise_post[v - 1]->InputFromMVN(resultMVNs.at(v - 1)->GetSubmatrix(
                m_num_params + 1, m_num_params + m_noise_params));
        }
        else
        {
            // Set the initial posterior for model params. Model
            // may want the voxel data in order to do this
            PassModelData(v);
            m_model->GetInitialPosterior(fwd_post[v - 1]);
            // Set initial noise posterior
            noise_post[v - 1] = initialNoisePosterior->Clone();
        }

        if (m_locked_linear)
        {
            lockedLinearCentres.Column(v)
                = lockedLinearDists.at(v - 1)->means.Rows(1, m_num_params);
            m_lin_model[v - 1].ReCentre(lockedLinearCentres.Column(v));
        }
        else
        {
            m_lin_model[v - 1].ReCentre(fwd_post[v - 1].means);
        }

        // Create per-voxel convergence detector. Initialization of m_needF is
        // inefficient but not harmful because all convergence detectors are the same type
        m_conv[v - 1] = ConvergenceDetector::NewFromName(conv_name);
        m_conv[v - 1]->Initialize(rundata);

        // Free energy
        m_saveF = rundata.GetBool("save-free-energy");
        m_printF = rundata.GetBool("print-free-energy");
        m_needF = m_conv[v - 1]->UseF() || m_printF || m_saveF;

        noise_prior[v - 1] = initialNoisePrior->Clone();
        m_noise->Precalculate(*noise_post[v - 1], *noise_prior[v - 1], m_origdata.Column(v));
    }
}

ThreadContext::~ThreadContext()
{
    while (!resultMVNs.empty())
    {
        delete resultMVNs.back();
        resultMVNs.pop_back();
    }
} 

void ThreadContext::IgnoreVoxel(int v)
{
    LOG << "Vb::IgnoreVoxel This voxel will be ignored in further updates" << endl;

    ignore_voxels.push_back(v);

    // Remove voxel from lists of neighbours of other voxels.
    // We identify affected voxels by looking in the neighbour
    // lists for the bad voxel, because any voxel which has
    // the bad voxel as a neighbour will be a neighbour of the
    // bad voxel
    vector<int> nn = neighbours[v - 1];
    for (vector<int>::iterator i = nn.begin(); i != nn.end(); ++i)
    {
        // Reference to list of neighbours of some other voxel which
        // has the bad voxel as a neighbour
        vector<int> &n2 = neighbours[*i - 1];

        n2.erase(std::remove(n2.begin(), n2.end(), v), n2.end());
    }

    // Same for next-nearest-neighbours
    nn = neighbours2[v - 1];
    for (vector<int>::iterator i = nn.begin(); i != nn.end(); ++i)
    {
        // Reference to list of neighbours of some other voxel which
        // has the bad voxel as a neighbour
        vector<int> &n2 = neighbours2[*i - 1];

        n2.erase(std::remove(n2.begin(), n2.end(), v), n2.end());
    }
}

void ThreadContext::CheckCoordMatrixCorrectlyOrdered()
{
    // Only 3D
    assert(m_coords.Nrows() == 3);

    // Go through each voxel one at a time apart from last
    for (int v = 1; v <= nvoxels - 1; v++)
    {
        // Find difference between current coords and next
        ColumnVector diff = m_coords.Column(v + 1) - m_coords.Column(v);

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
void ThreadContext::CalcNeighbours(int spatial_dims)
{
    if (nvoxels == 0)
        return;

    // Voxels must be ordered by increasing z, y and x values respectively
    // otherwise binary search for voxel by offset will not work
    CheckCoordMatrixCorrectlyOrdered();

    // Create a column vector with one entry per voxel.
    ColumnVector offsets(nvoxels);

    // Populate offsets with the offset into the
    // matrix of each voxel. We assume that co-ordinates
    // could be zero but not negative
    int xsize = m_coords.Row(1).Maximum() + 1;
    int ysize = m_coords.Row(2).Maximum() + 1;
    for (int v = 1; v <= nvoxels; v++)
    {
        int x = m_coords(1, v);
        int y = m_coords(2, v);
        int z = m_coords(3, v);
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
