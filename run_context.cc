
RunContext::RunContext(FwdModel *fwd_model, FabberRunData &rundata)
    : Loggable(rundata.GetLogger()), it(0), v(1), model(fwd_model), 
{
    
}

void RunContext::SetupPerVoxelDists(FabberRunData &rundata)
{
    // Initialized in voxel loop below (from file or default as required)
    noise_post.resize(m_nvoxels, NULL);
    noise_prior.resize(m_nvoxels, NULL);
    fwd_post.resize(m_nvoxels);

    // Re-centred in voxel loop below
    lin_model.resize(m_nvoxels, LinearizedFwdModel(model));

    // Static initialization for all voxels currently
    fwd_prior.resize(m_nvoxels, MVNDist(m_num_params, m_log));

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
