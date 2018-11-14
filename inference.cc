/*  inference.cc - General inference technique base class

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-20015 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference.h"
#include "easylog.h"

#include <newmatio.h>

#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace NEWMAT;

std::vector<std::string> InferenceTechnique::GetKnown()
{
    InferenceTechniqueFactory *factory = InferenceTechniqueFactory::GetInstance();
    return factory->GetNames();
}

InferenceTechnique *InferenceTechnique::NewFromName(const string &name)
{
    InferenceTechniqueFactory *factory = InferenceTechniqueFactory::GetInstance();
    InferenceTechnique *inf = factory->Create(name);
    if (inf == NULL)
    {
        throw InvalidOptionValue("method", name, "Unrecognized inference method");
    }
    return inf;
}

void InferenceTechnique::UsageFromName(const string &name, std::ostream &stream)
{
    stream << "Usage information for method: " << name << endl << endl;

    std::auto_ptr<InferenceTechnique> method(NewFromName(name));
    stream << method->GetDescription() << endl << endl << "Options: " << endl << endl;
    vector<OptionSpec> options;
    method->GetOptions(options);
    if (options.size() > 0)
    {
        for (vector<OptionSpec>::iterator iter = options.begin(); iter != options.end(); ++iter)
        {
            stream << *iter << endl;
        }
    }
}

InferenceTechnique::InferenceTechnique()
    : m_halt_bad_voxel(true)
{
}

void InferenceTechnique::Initialize(FabberRunData &rundata)
{
    m_log = rundata.GetLogger();
    m_debug = rundata.GetBool("debug");
    if (m_debug)
        LOG << setprecision(17);

    // Allow calculation to continue even with bad voxels
    // Note that this is a bad idea when spatialDims>0, because the bad voxel
    // will drag its neighbours around... but can never recover!  Maybe a more
    // sensible approach is to reset bad voxels to the prior on each iteration.
    // Spatial VB ignores this parameter, presumably for the above reason
    m_halt_bad_voxel = !rundata.GetBool("allow-bad-voxels");
    if (m_halt_bad_voxel)
    {
        LOG << "InferenceTechnique::Note: numerical errors in voxels will cause the program to "
               "halt.\n"
            << "InferenceTechnique::Use --allow-bad-voxels (with caution!) to keep on "
               "calculating.\n";
    }
    else
    {
        LOG << "InferenceTechnique::Using --allow-bad-voxels: numerical errors in a voxel will\n"
            << "InferenceTechnique::simply stop the calculation of that voxel.\n"
            << "InferenceTechnique::Check log for 'Going on to the next voxel' messages.\n"
            << "InferenceTechnique::Note that you should get very few (if any) exceptions like "
               "this;"
            << "InferenceTechnique::they are probably due to bugs or a numerically unstable model.";
    }
}

void InferenceTechnique::SaveResults(FabberRunData &rundata) const
{
    LOG << "InferenceTechnique::Preparing to save results..." << endl;
    int total_nvoxels = 0;
    for (int i=0; i<m_ctxs.size(); i++) {
        total_nvoxels += m_ctxs[i]->nvoxels;
    }

    vector<Parameter> params;
    m_ctxs[0]->m_model->GetParameters(rundata, params);
    int nparams = params.size();

    // Save the finalMVN NIFTI file
    if (rundata.GetBool("save-mvn"))
    {
        ColumnVector aOne(1);
        aOne = 1.0;
        Matrix mvns;

        for (int c=0; c<m_ctxs.size(); c++) {
            ThreadContext *ctx = m_ctxs[c];
            if (mvns.Nrows() == 0) {
                int mvn_nparams = nparams + ctx->m_noise_params;
                //cerr << mvn_nparams << ", " << nparams << ", " << ctx->m_noise_params << endl;
                mvns.ReSize(mvn_nparams * (mvn_nparams + 1) / 2 + mvn_nparams + 1, total_nvoxels);
            }
            for (int vox = 1; vox <= ctx->nvoxels; vox++)
            {
                MVNDist *mvn = ctx->resultMVNs[vox - 1];
                //cerr << mvn << endl;
                //cerr << mvn->means.Nrows() << ", " << mvns.Nrows() << endl;
                
                //cerr << nparams << ", " << vox << ", " << ctx->start_voxel + vox - 1 << ", " << mvn << endl;
                //cerr << mvns.Nrows() << ", " << mvns.Ncols() << endl; 
                //cerr << mvn->GetCovariance().AsColumn();
                //cerr << mvn->means.Nrows();
                mvns.Column(ctx->start_voxel + vox - 1)
                    = mvn->GetCovariance().AsColumn() & mvn->means & aOne;
            }
        }

        rundata.SaveVoxelData("finalMVN", mvns, VDT_MVN);
    }

    // Create individual files for each parameter's mean and Z-stat
    if (rundata.GetBool("save-mean") | rundata.GetBool("save-std") | rundata.GetBool("save-zstat"))
    {
        LOG << "InferenceTechnique::Writing means..." << endl;
        for (unsigned i = 1; i <= params.size(); i++)
        {
            Matrix paramMean, paramZstat, paramStd;
            paramMean.ReSize(1, total_nvoxels);
            paramZstat.ReSize(1, total_nvoxels);
            paramStd.ReSize(1, total_nvoxels);

            for (int c=0; c<m_ctxs.size(); c++) {
                ThreadContext *ctx = m_ctxs[c];
                for (int vox = 1; vox <= ctx->nvoxels; vox++)
                {
                    MVNDist *result = ctx->resultMVNs[vox - 1];
                    ctx->m_model->ToModel(*result);
                    paramMean(1, ctx->start_voxel + vox - 1) = result->means(i);
                    double std = sqrt(result->GetCovariance()(i, i));
                    paramZstat(1, ctx->start_voxel + vox - 1) = paramMean(1, ctx->start_voxel + vox - 1) / std;
                    paramStd(1, ctx->start_voxel + vox - 1) = std;
                }
            }

            if (rundata.GetBool("save-mean"))
                rundata.SaveVoxelData("mean_" + params.at(i - 1).name, paramMean);
            if (rundata.GetBool("save-zstat"))
                rundata.SaveVoxelData("zstat_" + params.at(i - 1).name, paramZstat);
            if (rundata.GetBool("save-std"))
                rundata.SaveVoxelData("std_" + params.at(i - 1).name, paramStd);
        }
    }

    // Produce the model fit and residual volume series
    bool saveModelFit = rundata.GetBool("save-model-fit");
    bool saveResiduals = rundata.GetBool("save-residuals");
    vector<string> outputs;
    outputs.push_back("");
    m_ctxs[0]->m_model->GetOutputs(outputs);
    if (saveModelFit || saveResiduals || (outputs.size() > 1))
    {
        LOG << "InferenceTechnique::Writing model time series data (fit, residuals and "
               "model-specific output)"
            << endl;

        Matrix result, residuals, datamtx, coords, suppdata;
        datamtx = rundata.GetMainVoxelData(); // it is just possible that the model needs the data
                                              // in its calculations
        coords = rundata.GetVoxelCoords();
        suppdata = rundata.GetVoxelSuppData();
        result.ReSize(datamtx.Nrows(), total_nvoxels);
        ColumnVector tmp;
        for (vector<string>::iterator iter = outputs.begin(); iter != outputs.end(); ++iter)
        {
            LOG << "InferenceTechnique::Evaluating model for output: " << *iter << endl;
            for (int i=0; i<m_ctxs.size(); i++) {
                ThreadContext *ctx = m_ctxs[i];
                for (int vox = 1; vox <= ctx->nvoxels; vox++)
                {
                    // do the evaluation
                    try
                    {
                        // pass in stuff that the model might need
                        ColumnVector y = datamtx.Column(vox);
                        ColumnVector vcoords = coords.Column(vox);
                        if (suppdata.Ncols() > 0)
                        {
                            ctx->m_model->PassData(vox, y, vcoords, suppdata.Column(vox));
                        }
                        else
                        {
                            ctx->m_model->PassData(vox, y, vcoords);
                        }

                        ctx->m_model->EvaluateFabber(
                            ctx->resultMVNs.at(vox - 1)->means.Rows(1, ctx->m_num_params), tmp, *iter);
                        if (result.Nrows() != tmp.Nrows())
                        {
                            // Only occurs on first voxel if output size is not equal to
                            // data size
                            result.ReSize(tmp.Nrows(), total_nvoxels);
                        }
                        result.Column(ctx->start_voxel + vox - 1) = tmp;
                    }
                    // Ignore exceptions for the default Evaluate key - errors when evaluating the model would already have
                    // occurred during inference and the relevant warnings output
                    catch (NEWMAT::Exception &e)
                    {
                        if (*iter != "") 
                        {
                            LOG << "InferenceTechnique::NEWMAT error generating output " << *iter << " for voxel " << vox << " : " << e.what() << endl;
                        } 
                    }
                    catch (Exception &e)
                    {
                        if (*iter != "") 
                        {
                            LOG << "InferenceTechnique::Error generating output " << *iter << " for voxel " << vox << " : " << e.what() << endl;
                        } 
                    }
                    catch (...)
                    {
                        if (*iter != "") 
                        {
                            LOG << "InferenceTechnique::Unexpected error generating output " << *iter << " for voxel " << vox << " : no message available" << endl;
                        }  
                    }
                }
            }
            if (*iter == "")
            {           
                if (saveResiduals)
                {
                    LOG << "InferenceTechnique::Saving residuals" << endl;
                    residuals = datamtx - result;
                    rundata.SaveVoxelData("residuals", residuals);
                }
                if (saveModelFit)
                {
                    LOG << "InferenceTechnique::Saving model fit" << endl;
                    rundata.SaveVoxelData("modelfit", result);
                }
            }
            else if (rundata.GetBool("save-model-extras"))
            {
                LOG << "InferenceTechnique::Saving extra output (size=" << result.Nrows() << ")" << endl;
                rundata.SaveVoxelData(*iter, result);
            }
        }
    }

#if 0
	{
		LOG << "InferenceTechnique::Writing model variances..." << endl;
		Matrix modelStd;
		Matrix datamtx = data.GetMainVoxelData();
		modelStd.ReSize(datamtx.Nrows(), nvoxels);
		for (int vox = 1; vox <= nvoxels; vox++)
		{
			LinearizedFwdModel lin(model);
			lin.ReCentre(resultMVNs.at(vox - 1)->means.Rows(1, m_num_params));
			Matrix var = resultMVNs.at(vox - 1)->GetCovariance().SymSubMatrix(1, m_num_params);

			Matrix mvar = lin.Jacobian() * var * lin.Jacobian().t();
			ColumnVector tmp(datamtx.Nrows());
			for (int i = 1; i <= datamtx.Nrows(); i++)
			{
				tmp(i) = sqrt(mvar(i, i));
			}

			modelStd.Column(vox) = tmp;
		}
		data.SaveVoxelData("modelstd", modelStd);

	}
#endif

    LOG << "InferenceTechnique::Done writing results." << endl;
}

InferenceTechnique::~InferenceTechnique()
{
    //if (m_ctx) {
    //    delete m_ctx;
   // } FIXME
}
