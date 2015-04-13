/*  inference.h - General inference technique base class

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */


#if !defined(__inference_h)
#define __inference_h


#pragma once
#include <map>
#include <string>
#include <vector>
#include "fwdmodel.h"
#include "noisemodel.h"
#include "easylog.h"
#include "easyoptions.h"
#include "dataset.h"
#ifdef __FABBER_MOTION
 #include "Update_deformation.h"
 #include "mcflirt/rigidreglib.h"
#endif //__FABBER_MOTION

class InferenceTechnique {
    
 public:
  static InferenceTechnique* NewFromName(const string& method);
    // returns a *new* instance of an Inference-Technique-derived class,
    // as determined by the name given in "method".
    
 public:
  InferenceTechnique() : model(NULL), noise(NULL) { return; }
  virtual void Setup(ArgsType& args);
  virtual void SetOutputFilenames(const string& output)
    { outputDir = output; }
  virtual void DoCalculations(const DataSet& data) = 0;
  virtual void SaveResults(const DataSet& data) const;
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

#endif

