/* setup.cc - FABBER default configuration setup class declarations.

   Mike Jackson, The University of Edinburgh & Michael Chappell, FMRIB Analysis Group & IBME QuBIc group

   Copyright (C) 2015 University of Oxford  */


#include "inference_nlls.h"
#include "inference_spatialvb.h"
#include "inference_vb.h"

#include "noisemodel_ar.h"
#include "noisemodel_white.h"

#include "setup.h"

void FabberSetup::SetupDefaultInferenceTechniques() {
  InferenceTechniqueFactory* factory = 
    InferenceTechniqueFactory::GetInstance();
  factory->Add("vb", &VariationalBayesInferenceTechnique::NewInstance);
  factory->Add("spatialvb", &SpatialVariationalBayes::NewInstance);
  factory->Add("nlls", &NLLSInferenceTechnique::NewInstance);
}

void FabberSetup::SetupDefaultNoiseModels() {
  NoiseModelFactory* factory = NoiseModelFactory::GetInstance();
  factory->Add("ar", &Ar1cNoiseModel::NewInstance);
  factory->Add("white", &WhiteNoiseModel::NewInstance);
}

void FabberSetup::SetupDefaultFwdModels() {
  FwdModelFactory* factory = FwdModelFactory::GetInstance();
  factory->Add("linear", &LinearFwdModel::NewInstance);
}

void FabberSetup::SetupDefaults() {
  FabberSetup::SetupDefaultInferenceTechniques();
  FabberSetup::SetupDefaultNoiseModels();
  FabberSetup::SetupDefaultFwdModels();
}

void FabberSetup::Destroy() {
  FwdModelFactory::Destroy();
  NoiseModelFactory::Destroy();
  InferenceTechniqueFactory::Destroy();
}
