// fwdmodel_sine.h - A simple sine curve fitting model

#include "fwdmodel.h"
#include "dataset.h"
#include <string>
#include "newmat.h"

class SineFwdModel: public FwdModel
{
public:
	static FwdModel* NewInstance();
	void GetOptions(std::vector<OptionSpec> &opts) const;

	void Initialize(FabberRunData& args);
	int NumParams() const;
	void NameParams(std::vector<std::string>& names) const;
	string ModelVersion() const;
	void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;
	void Evaluate(const NEWMAT::ColumnVector& params, NEWMAT::ColumnVector& result) const;
private:
	bool m_include_offset;
	static FactoryRegistration<FwdModelFactory, SineFwdModel> registration;
};

