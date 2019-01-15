//  fwdmodel_sine.cc - Implements a simple sine curve fitting model
#include "fwdmodel_sine.h"

#include "fabber_core/fwdmodel.h"

#include <math.h>

using namespace std;
using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, SineFwdModel> SineFwdModel::registration("sine");

FwdModel *SineFwdModel::NewInstance()
{
    return new SineFwdModel();
}

string SineFwdModel::ModelVersion() const
{
    return "1.0";
}

string SineFwdModel::GetDescription() const
{
    return "Example model which uses a sine function";
}

static OptionSpec OPTIONS[] = {
    { "use-offset", OPT_BOOL, "If True, allow an additional constant offset parameter", OPT_NONREQ, "false" },
    { "" }
};

void SineFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

void SineFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    int p=0;
    params.push_back(Parameter(p++, "a", DistParams(1, 1e6), DistParams(1, 1e6)));
    params.push_back(Parameter(p++, "b", DistParams(1, 1e6), DistParams(1, 1e6)));
    params.push_back(Parameter(p++, "c", DistParams(0, 1e6), DistParams(0, 1e6)));
    if (m_include_offset) {
        params.push_back(Parameter(p++, "d", DistParams(0, 1e6), DistParams(0, 1e6)));
    }
}

void SineFwdModel::Initialize(FabberRunData &rundata)
{
    m_include_offset = rundata.GetBool("use-offset");
}

void SineFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Check we have been given the right number of parameters
    assert(params.Nrows() == NumParams());
    float a = params(1);
    float b = params(2);
    float c = params(3);
    result.ReSize(data.Nrows());

    for (int i = 1; i <= data.Nrows(); i++)
    {
        float t = float(i) / data.Nrows();
        double res = a * sin(b * (t - c));
        if (m_include_offset)
            res += params(4);
        result(i) = res;
    }
}
