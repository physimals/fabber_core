//  fwdmodel_exp.cc - Implements a simple exponential decay fitting model
#include "fwdmodel_exp.h"

#include <fabber_core/fwdmodel.h>
#include <fabber_core/priors.h>

#include <math.h>

using namespace std;
using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, ExpFwdModel> ExpFwdModel::registration("exp");

FwdModel *ExpFwdModel::NewInstance()
{
    return new ExpFwdModel();
}

string ExpFwdModel::ModelVersion() const
{
    return "1.0";
}

string ExpFwdModel::GetDescription() const
{
    return "Example model of a sum of exponentials";
}

static OptionSpec OPTIONS[] = {
    { "dt", OPT_FLOAT, "Time separation between samples", OPT_REQ, "" },
    { "num-exps", OPT_INT, "Number of independent decay rates", OPT_NONREQ, "1" },
    { "" }
};

void ExpFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

void ExpFwdModel::Initialize(FabberRunData& rundata)
{
    m_dt = rundata.GetDouble("dt");
    m_num = rundata.GetIntDefault("num-exps", 1);
}

void ExpFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    int p=0;
    for (int i=0; i<m_num; i++) {
        params.push_back(Parameter(p++, "amp" + stringify(i+1), DistParams(1, 100), DistParams(1, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
        params.push_back(Parameter(p++, "r" + stringify(i+1), DistParams(1, 100), DistParams(1, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
        //params.push_back(Parameter(p++, "amp" + stringify(i+1), DistParams(1.0, 1e6), DistParams(1, 100)));
        //params.push_back(Parameter(p++, "r" + stringify(i+1), DistParams(1, 1e6), DistParams(1, 100)));
    }
}

void ExpFwdModel::EvaluateModel(const NEWMAT::ColumnVector &params, 
                                NEWMAT::ColumnVector &result, 
                                const std::string &key) const
{
    result.ReSize(data.Nrows());
    result = 0;
    
    for (int i=0; i<m_num; i++) {
        double amp = params(2*i+1);
        double r = params(2*i+2);
        //if (amp < 1e-6) amp = 1e-6;
        //if (r < 1e-6) r = 1e-6;
        //if (amp > 100) amp = 100;
        //if (r > 100) r = 100;
        for (int i=0; i < data.Nrows(); i++)
        {
            double t = double(i) * m_dt;
            double val = amp * exp(-r * t);
            result(i+1) += val;
        }
    }
}

void ExpFwdModel::InitVoxelPosterior(MVNDist &posterior) const
{
    double data_max = data.Maximum();

    for (int i=0; i<m_num; i++) {
        posterior.means(2*i+1) = data_max / (m_num+i);
    }
}
