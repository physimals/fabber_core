//  fwdmodel_sine.cc - Implements a simple sine curve fitting model

#include "fwdmodel_sine.h"

#include "fabber_core/fwdmodel.h"

#include <math.h>

using namespace std;
using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, SineFwdModel> SineFwdModel::registration("sine");

FwdModel* SineFwdModel::NewInstance()
{
    return new SineFwdModel();
}

string SineFwdModel::GetDescription() const
{
    return "Example model which uses a sine function";
}

string SineFwdModel::ModelVersion() const
{
    return "1.0";
}

static OptionSpec OPTIONS[] = {
    { "use-offset", OPT_BOOL, "If True, allow an additional constant offset parameter", OPT_NONREQ, "false" },
    { "" }
};

void SineFwdModel::GetOptions(vector<OptionSpec>& opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++) {
        opts.push_back(OPTIONS[i]);
    }
}

void SineFwdModel::Initialize(FabberRunData& rundata)
{
    m_include_offset = rundata.GetBool("use-offset");
}

int SineFwdModel::NumParams() const
{
    if (m_include_offset)
        return 4;
    else
        return 3;
}

void SineFwdModel::NameParams(vector<string>& names) const
{
    names.clear();
    names.push_back("a");
    names.push_back("b");
    names.push_back("c");
    if (m_include_offset)
        names.push_back("d");
}

void SineFwdModel::HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const
{
    int num_params = NumParams();
    // Check we have been given a distribution of the right number of parameters
    assert(prior.means.Nrows() == num_params);
    prior.means = 0;
    prior.SetPrecisions(IdentityMatrix(num_params) * 1e-12);
    posterior = prior;
}

void SineFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
    // Check we have been given the right number of parameters
    assert(params.Nrows() == NumParams());
    result.ReSize(data.Nrows());

    for (int i = 1; i <= data.Nrows(); i++) {
        double res = params(1) * sin(params(2) * (i - params(3)));
        if (m_include_offset)
            res += params(4);
        result(i) = res;
    }
}
