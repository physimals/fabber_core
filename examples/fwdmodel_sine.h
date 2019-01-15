// fwdmodel_sine.h - A simple sine curve fitting model
#pragma once

#include "fabber_core/fwdmodel.h"

#include "newmat.h"

#include <string>
#include <vector>

class SineFwdModel : public FwdModel {
public:
    static FwdModel* NewInstance();

    SineFwdModel()
        : m_include_offset(false)
    {
    }

    std::string ModelVersion() const;
    std::string GetDescription() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;

    void Initialize(FabberRunData &args);
    void EvaluateModel(const NEWMAT::ColumnVector &params, 
                       NEWMAT::ColumnVector &result, 
                       const std::string &key="") const;
    
private:
    bool m_include_offset;
    static FactoryRegistration<FwdModelFactory, SineFwdModel> registration;
};

