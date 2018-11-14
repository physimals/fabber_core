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

void *ThreadContext::Start(void *obj)
{
    ThreadContext *tc  = ((ThreadContext *)obj);
    tc->Run();
    return NULL;
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
/* FIXME
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
    }*/
}
