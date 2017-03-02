// fwdmodel_sine.h - A simple sine curve fitting model

#ifndef FWDMODEL_SINE_H
#define FWDMODEL_SINE_H

#include "fabber_core/fwdmodel.h"

#include "newmat.h"

#include <string>

class SineFwdModel : public FwdModel {
public:
    static FwdModel* NewInstance();

    SineFwdModel()
        : m_include_offset(false)
    {
    }

    void GetOptions(std::vector<OptionSpec>& opts) const;
    std::string GetDescription() const;
    std::string ModelVersion() const;

    void Initialize(FabberRunData& args);
    int NumParams() const;
    void NameParams(std::vector<std::string>& names) const;
    void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;
    void Evaluate(const NEWMAT::ColumnVector& params, NEWMAT::ColumnVector& result) const;

private:
    bool m_include_offset;
    static FactoryRegistration<FwdModelFactory, SineFwdModel> registration;
};

#endif
