/*  inference.h - General inference technique base class

    Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

    Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */


#ifndef __FABBER_INFERENCE_H
#define __FABBER_INFERENCE_H 1

#include <map>
#include <string>
#include <typeinfo>
#include <vector>

#include "dataset.h"
#include "fwdmodel.h"
#include "easylog.h"
#include "easyoptions.h"
#include "noisemodel.h"
#include "utils.h"

#ifdef __FABBER_MOTION
 #include "Update_deformation.h"
 #include "mcflirt/rigidreglib.h"
#endif //__FABBER_MOTION

class InferenceTechnique {
        
 public:
  /**
     * Create a new instance of this class.
     * @return pointer to new instance.
     */
    static InferenceTechnique* NewInstance();
    /**
     * Constructor.
     */
    InferenceTechnique() : model(NULL), noise(NULL) { 
      // No-op.
    }
    /**
     * Initialize a new instance to use the given forward model
     * and extract additional configuration from the given
     * arguments.
     * @param fwd_model Forward model to be used.
     * @param args Additional configuration parameters.
     */
    virtual void Initialize(FwdModel* fwd_model, ArgsType& args);
    /**
     * Set the path to the output directory.
     * @param output Directory path, can be relative or 
     * absolute.
     */
    virtual void SetOutputFilenames(const string& output) { 
      outputDir = output; 
    }
    /**
     * Perform inference using the given model upon the given data.
     * @param data 
     */
    virtual void DoCalculations(const DataSet& data) = 0;
    /**
     * Save the results in the output directory.
     */
    virtual void SaveResults(const DataSet& data) const;
    /**
     * Destructor.
     */
    virtual ~InferenceTechnique();

 protected:
  FwdModel* model;
  NoiseModel* noise;
  string outputDir;
  bool saveModelFit;
  bool saveResiduals;
  
  vector<MVNDist*> resultMVNs;
  vector<MVNDist*> resultMVNsWithoutPrior; // optional; used by Adrian's spatial priors research
  vector<double> resultFs;

  void InitMVNFromFile(vector<MVNDist*>& continueFromDists,string continueFromFile, const DataSet& allData, string paramFilename);
  
  // Motion related stuff
  int Nmcstep; // number of motion correction steps to run

private:
    const InferenceTechnique& operator=(const InferenceTechnique& from)
        { assert(false); return from; } // just not allowed. 

};

/** 
 * \ref SingletonFactory that returns pointers to 
 * \ref InferenceTechnique.
 */
typedef SingletonFactory<InferenceTechnique> InferenceTechniqueFactory;

// Motion Correction class
#ifdef __FABBER_MOTION
//   NB: for now the mask should cover the *entire* image as we zero everything
//       outside of the mask, which is not good for registration
//       In future we'd need allData to be able to provide the original image (or something to)

class MCobj {
public:
  MCobj(const DataSet& allData, int dof);
  void run_mc(const Matrix& modelpred_mat, Matrix& finalimage_mat);
  void set_num_iter(int nit) { num_iter=nit; }
private:
  int userdof;  // anything over 13 is full nonlinear
  int num_iter;  // default 10
  volume<float> mask;
  Matrix affmat;
  mcflirt mcf;
  volume4D<float> defx;
  volume4D<float> defy;
  volume4D<float> defz;
  // things below are kept for efficiency (?) in order to avoid repeated allocation/destruction
  volume4D<float> tmpx;
  volume4D<float> tmpy;
  volume4D<float> tmpz;
  volume4D<float> modelpred;
  volume4D<float> finalimage;
  volume4D<float> wholeimage;
};

#endif // __FABBER_MOTION

#endif // __FABBER_INFERENCE_H

