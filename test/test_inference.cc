//
// Test for general inference methods
//
// Separate classes are used to test specific functionality of
// particular inference methods

#include "gtest/gtest.h"

#include "inference.h"
#include "dataset.h"
#include "setup.h"
#include "easylog.h"

namespace {

// The fixture for testing class Foo.
class InferenceMethodTest : public ::testing::TestWithParam<string> {
 protected:
  InferenceMethodTest() {
    FabberSetup::SetupDefaults();
    EasyLog::StartLog(".", true);
  }

  virtual ~InferenceMethodTest() {
    FabberSetup::Destroy();
     EasyLog::StopLog();
 }

  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};

// Tests that the VB inference method can be created
TEST_P(InferenceMethodTest, CanCreate) {
  InferenceTechnique *vb = InferenceTechnique::NewFromName(GetParam());
  ASSERT_TRUE(vb);
}

// Tests that the VB inference method on a trivial one-parameter model
// with 1 voxel, 1 timeslice
TEST_P(InferenceMethodTest, DISABLED_OneParamOneVoxelOneTimeslice) 
{
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data;
    data.ReSize(1, 1);
    voxelCoords.ReSize(3, 1);
    data(1, 1) = VAL;
    voxelCoords << 1 << 1 << 1;

    FabberRunData rundata;
    rundata.SetVoxelCoords(voxelCoords);
    rundata.SetMainVoxelData(data);
    rundata.Set("noise", "white");
    rundata.Set("model", "trivial");
    rundata.SetBool("print-free-energy", true);
    rundata.Set("method", GetParam());
    ASSERT_NO_THROW(rundata.Run());

    NEWMAT::Matrix mean = rundata.GetVoxelData("mean_p");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), 1);
    ASSERT_FLOAT_EQ(mean(1, 1), VAL);
}

// Tests that the VB inference method on a trivial one-parameter model
// with 1 voxel and multiple timeslices
TEST_P(InferenceMethodTest, OneParamOneVoxelMultiTimeslice) 
{
    int NTIMES = 10;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data;
    data.ReSize(NTIMES, 1);
    voxelCoords.ReSize(3, 1);
    voxelCoords << 1 << 1 << 1;
    for (int i=0; i<NTIMES; i++) {
      data(i+1, 1) = VAL;
    }

    FabberRunData rundata;
    rundata.SetVoxelCoords(voxelCoords);
    rundata.SetMainVoxelData(data);
    rundata.SetBool("print-free-energy", true);
    rundata.Set("noise", "white");
    rundata.Set("model", "trivial");
    rundata.Set("method", GetParam());
    rundata.Run();

    NEWMAT::Matrix mean = rundata.GetVoxelData("mean_p");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), 1);
    ASSERT_FLOAT_EQ(mean(1, 1), VAL);
}

// Tests that the VB inference method on a trivial one-parameter model
// with multiple voxels and multiple timeslices
TEST_P(InferenceMethodTest, OneParamMultiVoxelMultiTimeslice) 
{
    int NTIMES = 10;
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data;
    data.ReSize(NTIMES, VSIZE*VSIZE*VSIZE);
    voxelCoords.ReSize(3, VSIZE*VSIZE*VSIZE);
    int v=1;
    for (int z = 0; z < VSIZE; z++) {
        for (int y = 0; y < VSIZE; y++) {
            for (int x = 0; x < VSIZE; x++) {
                voxelCoords(1, v) = x;
                voxelCoords(2, v) = y;
                voxelCoords(3, v) = z;
                for (int n = 0; n < NTIMES; n++) {
                    data(n + 1, v) = VAL;
                }
                v++;
            }
        }
    }

    FabberRunData rundata;
    rundata.SetVoxelCoords(voxelCoords);
    rundata.SetMainVoxelData(data);
    rundata.Set("noise", "white");
    rundata.Set("model", "trivial");
    rundata.Set("method", GetParam());
    rundata.Run();

    NEWMAT::Matrix mean = rundata.GetVoxelData("mean_p");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE*VSIZE*VSIZE);
    for (int i=0; i<VSIZE*VSIZE*VSIZE; i++) {
      ASSERT_FLOAT_EQ(mean(1, i+1), VAL);
  }
}

// Tests that the VB inference method on a trivial one-parameter model
// with multiple voxels and multiple timeslices
TEST_P(InferenceMethodTest, OneParamMultiVoxelMultiTimesliceVariable) 
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data;
    data.ReSize(NTIMES, VSIZE*VSIZE*VSIZE);
    voxelCoords.ReSize(3, VSIZE*VSIZE*VSIZE);
    int v=1;
    for (int z = 0; z < VSIZE; z++) {
        for (int y = 0; y < VSIZE; y++) {
            for (int x = 0; x < VSIZE; x++) {
                voxelCoords(1, v) = x;
                voxelCoords(2, v) = y;
                voxelCoords(3, v) = z;
                for (int n = 0; n < NTIMES; n++) {
                    if (n % 2 == 0) data(n + 1, v) = VAL;
                    else data(n + 1, v) = VAL*3;
                }
                v++;
            }
        }
    }

    FabberRunData rundata;
    rundata.SetVoxelCoords(voxelCoords);
    rundata.SetMainVoxelData(data);
    rundata.Set("noise", "white");
    rundata.Set("model", "trivial");
    rundata.Set("method", GetParam());
    rundata.Run();

    NEWMAT::Matrix mean = rundata.GetVoxelData("mean_p");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE*VSIZE*VSIZE);
    for (int i=0; i<VSIZE*VSIZE*VSIZE; i++) {
      ASSERT_FLOAT_EQ(mean(1, i+1), VAL*2);
  }
}

INSTANTIATE_TEST_CASE_P(MethodTests,
                        InferenceMethodTest,
                        ::testing::Values("vb", "nlls", "spatialvb"));


}  // namespace

