/*  inference_vb.h - VB inference technique class declarations

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

    Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once
#include "inference.h"

// Subclasses defined here!
// The actual implementation details don't need to be visible to anyone else,
// since the appropriate derived class is returned by the global function
// PickInferenceTechnique(string).

// Forward declaration -- see convergence.h, which is only actually
// #included in inference_vb.cc.  For now it's good enough just to know
// that the class exists.
class ConvergenceDetector;

class VariationalBayesInferenceTechnique : public InferenceTechnique {
   public:
  static InferenceTechnique* NewInstance();
      VariationalBayesInferenceTechnique() : conv(NULL), 
        initialFwdPrior(NULL), initialFwdPosterior(NULL), 
        initialNoisePrior(NULL), initialNoisePosterior(NULL) { return; }
      virtual void Initialize(FwdModel* model, ArgsType& args);
      //  virtual void SetOutputFilenames(ArgsType& args);
      virtual void DoCalculations(const DataSet& data);    
      virtual ~VariationalBayesInferenceTechnique();
   protected:
      ConvergenceDetector* conv;
      const MVNDist* initialFwdPrior;
      const MVNDist* initialFwdPosterior;
      const NoiseParams* initialNoisePrior;
      const NoiseParams* initialNoisePosterior;

      // specificaiont of priors from command line
      vector<unsigned int> PSPidx;
      string PriorsTypes;
      vector<string> imagepriorstr;
      
      // These are used for resuming a previous calculation
      string continueFromFile; // if empty, use initial posterior dists above
      string paramFilename;
      bool continueFwdOnly; // Only have fwd-model information
      
      // Reduce this to a linear problem, using the given 
      // voxelwise linearizations (probably loaded from an MVN)
      string lockedLinearFile; 

      bool haltOnBadVoxel;
      bool printF;
      bool needF;
};

