/* setup.cc - FABBER default configuration setup class declarations.

 Mike Jackson, The University of Edinburgh & Michael Chappell, FMRIB Analysis Group & IBME QuBIc
 group

 Copyright (C) 2015 University of Oxford  */

#include "setup.h"

#include "inference.h"
#include "inference_vb.h"
#ifndef NO_NLLS
#include "inference_nlls.h"
#endif

#include "noisemodel.h"
#include "noisemodel_ar.h"
#include "noisemodel_white.h"

#include "fwdmodel.h"
#include "fwdmodel_linear.h"
#include "fwdmodel_poly.h"

#include "convergence.h"

void FabberSetup::SetupDefaultInferenceTechniques()
{
    InferenceTechniqueFactory *factory = InferenceTechniqueFactory::GetInstance();
    factory->Add("vb", &Vb::NewInstance);
    factory->Add("spatialvb", &Vb::NewInstance);
#ifndef NO_NLLS
    factory->Add("nlls", &NLLSInferenceTechnique::NewInstance);
#endif
}

void FabberSetup::SetupDefaultNoiseModels()
{
    NoiseModelFactory *factory = NoiseModelFactory::GetInstance();
    factory->Add("ar", &Ar1cNoiseModel::NewInstance);
    factory->Add("white", &WhiteNoiseModel::NewInstance);
}

void FabberSetup::SetupDefaultFwdModels()
{
    FwdModelFactory *factory = FwdModelFactory::GetInstance();
    factory->Add("linear", &LinearFwdModel::NewInstance);
    factory->Add("poly", &PolynomialFwdModel::NewInstance);
}

void FabberSetup::SetupDefaultConvergenceDetectors()
{
    ConvergenceDetectorFactory *factory = ConvergenceDetectorFactory::GetInstance();
    factory->Add("maxits", &CountingConvergenceDetector::NewInstance);
    factory->Add("pointzeroone", &FchangeConvergenceDetector::NewInstance);
    factory->Add("freduce", &FreduceConvergenceDetector::NewInstance);
    factory->Add("trialmode", &TrialModeConvergenceDetector::NewInstance);
    factory->Add("lm", &LMConvergenceDetector::NewInstance);
}

void FabberSetup::SetupDefaults()
{
    FabberSetup::SetupDefaultInferenceTechniques();
    FabberSetup::SetupDefaultNoiseModels();
    FabberSetup::SetupDefaultFwdModels();
    FabberSetup::SetupDefaultConvergenceDetectors();
}

void FabberSetup::Destroy()
{
    FwdModelFactory::Destroy();
    NoiseModelFactory::Destroy();
    InferenceTechniqueFactory::Destroy();
}
