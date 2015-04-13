/*  dataset.h - Data-loading class for FABBER

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "easyoptions.h"
#include "easylog.h"
#include <vector>

#ifndef __FABBER_LIBRARYONLY
#include "newimage/newimageall.h"
#endif //__FABBER_LIBRARYONLY

using namespace Utilities;

class DataSet
{
 public:
  void LoadData(ArgsType& args);

#ifndef __FABBER_LIBRARYONLY
  const NEWIMAGE::volume<float>& GetMask() const { return mask; }
#endif // __FABBER_LIBRARYONLY
  const NEWMAT::Matrix& GetVoxelData() const { return voxelData; }
  const NEWMAT::Matrix& GetVoxelCoords() const { return voxelCoords; }
  const NEWMAT::Matrix& GetVoxelSuppData() const { return voxelSuppData; }

 protected:
#ifndef __FABBER_LIBRARYONLY
  NEWIMAGE::volume<float> mask; // Will be unset if UsingMatrixIO!
#endif //__FABBER_LIBRARYONLY
  NEWMAT::Matrix voxelData;

  // supplementary data (timeseries)
  NEWMAT::Matrix voxelSuppData;

  // coordinates of each voxel
  NEWMAT::Matrix voxelCoords;  // is 3 x Nvox; integer indices (from 0), NOT mm positions
};



