/* setup.h - FABBER default configuration setup class declarations.

 Mike Jackson, The University of Edinburgh & Michael Chappell, FMRIB Analysis Group & IBME QuBIc group

 Copyright (C) 2015 University of Oxford  */

/*  CCOPYRIGHT */
#pragma once

/**
 * Functions to record the available inference techniques, noise
 * models and forward models so they are available across FABBER.
 */
class FabberSetup {
public:
    /**
	 * Invoke each Setup function in turn.
	 */
    static void SetupDefaults();
    /**
	 * Populate singleton \ref InferenceTechniqueFactory with default
	 * inference techniques.
	 */
    static void SetupDefaultInferenceTechniques();
    /**
	 * Populate singleton \ref NoiseModelFactory with default noise
	 * models.
	 */
    static void SetupDefaultNoiseModels();
    /**
	 *
	 * Populate singleton \ref FwdModelFactory with default forward
	 * models.
	 */
    static void SetupDefaultFwdModels();

    /**
	 *
	 * Populate singleton \ref ConvergenceDetectorFactory with default
	 * convergence detectors
	 */
    static void SetupDefaultConvergenceDetectors();
    /**
	 * Destroy all singleton factory instances.
	 */
    static void Destroy();
};
