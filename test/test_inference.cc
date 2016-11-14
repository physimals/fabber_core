//
// Generic tests to apply to all inference methods
//
// Separate classes are used to test specific functionality of
// particular inference methods

#include "gtest/gtest.h"

#include "inference.h"
#include "dataset.h"
#include "setup.h"
#include "easylog.h"

#include <fstream>

namespace
{

// The fixture for testing class Foo.
class InferenceMethodTest: public ::testing::TestWithParam<string>
{
protected:
	InferenceMethodTest()
	{
		FabberSetup::SetupDefaults();
		EasyLog::StartLog(".", true);
	}

	virtual ~InferenceMethodTest()
	{
		FabberSetup::Destroy();
		EasyLog::StopLog();
	}

	virtual void SetUp()
	{
	}

	virtual void TearDown()
	{
	}

	bool FloatEq(double d1, double d2, double epsilon = 0.001)
	{
		double diff = d1 - d2;
		if (diff < 0)
			diff = -diff;
		return diff < epsilon;
	}
};

// Tests that the VB inference method can be created
TEST_P(InferenceMethodTest, CanCreate)
{
	InferenceTechnique *vb = InferenceTechnique::NewFromName(GetParam());
	ASSERT_TRUE(vb);
}

// Tests inference on a trivial one-parameter model
// with 1 voxel, 1 timeslice. This is just for NLLs, VB cannot
// do this because with 1 model param + noise the model is overspecified
// for the data.
TEST_F(InferenceMethodTest, OneParamOneVoxelOneTimeslice)
{
	float VAL = 7.32;

	// Create coordinates and data matrices
	NEWMAT::Matrix voxelCoords, data;
	data.ReSize(1, 1);
	voxelCoords.ReSize(3, 1);
	data(1, 1) = VAL;
	voxelCoords << 1 << 1 << 1;

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data", data);
	rundata.Set("noise", "white");
	rundata.Set("model", "trivial");
	rundata.SetBool("print-free-energy", true);
	rundata.Set("method", "nlls");
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
	for (int i = 0; i < NTIMES; i++)
	{
		data(i + 1, 1) = VAL;
	}

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data", data);
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
	data.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	voxelCoords.ReSize(3, VSIZE * VSIZE * VSIZE);
	int v = 1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				for (int n = 0; n < NTIMES; n++)
				{
					data(n + 1, v) = VAL;
				}
				v++;
			}
		}
	}

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data", data);
	rundata.Set("noise", "white");
	rundata.Set("model", "trivial");
	rundata.Set("method", GetParam());
	rundata.Run();

	NEWMAT::Matrix mean = rundata.GetVoxelData("mean_p");
	ASSERT_EQ(mean.Nrows(), 1);
	ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
	for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
	{
		ASSERT_FLOAT_EQ(mean(1, i + 1), VAL);
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
	data.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	voxelCoords.ReSize(3, VSIZE * VSIZE * VSIZE);
	int v = 1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				for (int n = 0; n < NTIMES; n++)
				{
					if (n % 2 == 0)
						data(n + 1, v) = VAL;
					else
						data(n + 1, v) = VAL * 3;
				}
				v++;
			}
		}
	}

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data", data);
	rundata.Set("noise", "white");
	rundata.Set("model", "trivial");
	rundata.Set("method", GetParam());
	rundata.Run();

	NEWMAT::Matrix mean = rundata.GetVoxelData("mean_p");
	ASSERT_EQ(mean.Nrows(), 1);
	ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
	for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
	{
		ASSERT_FLOAT_EQ(mean(1, i + 1), VAL * 2);
	}
}

// Tests that we can read config from a file
TEST_P(InferenceMethodTest, ConfigFile)
{
	int NTIMES = 10; // needs to be even
	int VSIZE = 5;
	float VAL = 7.32;
	string FILENAME = "test_config";

	// Create coordinates and data matrices
	NEWMAT::Matrix voxelCoords, data;
	data.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	voxelCoords.ReSize(3, VSIZE * VSIZE * VSIZE);
	int v = 1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				for (int n = 0; n < NTIMES; n++)
				{
					if (n % 2 == 0)
						data(n + 1, v) = VAL;
					else
						data(n + 1, v) = VAL * 3;
				}
				v++;
			}
		}
	}

	ofstream os;
	os.open(FILENAME.c_str(), ios::out);
	os << "--noise=white" << endl << "--model=trivial" << endl << "--method=" << GetParam() << endl;
	os.close();

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data", data);
	rundata.ParseOldStyleParamFile(FILENAME);
	rundata.Run();

	NEWMAT::Matrix mean = rundata.GetVoxelData("mean_p");
	ASSERT_EQ(mean.Nrows(), 1);
	ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
	for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
	{
		ASSERT_FLOAT_EQ(mean(1, i + 1), VAL * 2);
	}
	remove(FILENAME.c_str());
}

// Tests that we can read config from command line arguments
TEST_P(InferenceMethodTest, CLArgs)
{
	int NTIMES = 10; // needs to be even
	int VSIZE = 5;
	float VAL = 7.32;
	string FILENAME = "test_config";

	// Create coordinates and data matrices
	NEWMAT::Matrix voxelCoords, data;
	data.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
	voxelCoords.ReSize(3, VSIZE * VSIZE * VSIZE);
	int v = 1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				for (int n = 0; n < NTIMES; n++)
				{
					if (n % 2 == 0)
						data(n + 1, v) = VAL;
					else
						data(n + 1, v) = VAL * 3;
				}
				v++;
			}
		}
	}

	string method = (string) "--method=" + GetParam();
	const char *argv[] =
	{ "fabber", "--noise=white", "--model=trivial", method.c_str() };
	int argc = 4;

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data", data);
	rundata.Parse(argc, (char **) argv);
	ASSERT_NO_THROW(rundata.Run());

	NEWMAT::Matrix mean = rundata.GetVoxelData("mean_p");
	ASSERT_EQ(mean.Nrows(), 1);
	ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
	for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
	{
		ASSERT_FLOAT_EQ(mean(1, i + 1), VAL * 2);
	}
	remove(FILENAME.c_str());
}

// Test fitting to a simple polynomial model
TEST_P(InferenceMethodTest, PolynomialFit)
{
	int NTIMES = 10;
	int VSIZE = 5;
	float VAL = 2;
	int DEGREE = 3;
	int n_voxels = VSIZE * VSIZE * VSIZE;

	// Create coordinates and data matrices
	// Data fitted to a cubic function
	NEWMAT::Matrix voxelCoords, data;
	data.ReSize(NTIMES, n_voxels);
	voxelCoords.ReSize(3, n_voxels);
	int v = 1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				for (int n = 0; n < NTIMES; n++)
				{
					// function is VAL + 1.5VAL x n^2 - 2*VAL*n^3
					data(n + 1, v) = VAL + (1.5 * VAL) * (n + 1) * (n + 1) - 2 * VAL * (n + 1) * (n + 1) * (n + 1);
				}
				v++;
			}
		}
	}

	// Do just 1 iteration

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data", data);
	rundata.Set("noise", "white");
	rundata.Set("model", "poly");
	rundata.Set("max-iterations", "50");
	rundata.Set("degree", stringify(DEGREE));
	rundata.Set("method", GetParam());
	rundata.Run();

	NEWMAT::Matrix mean = rundata.GetVoxelData("mean_c0");
	ASSERT_EQ(mean.Nrows(), 1);
	ASSERT_EQ(mean.Ncols(), n_voxels);
	for (int i = 0; i < n_voxels; i++)
	{
		ASSERT_TRUE(FloatEq(VAL, mean(1, i + 1)));
	}

	mean = rundata.GetVoxelData("mean_c1");
	ASSERT_EQ(mean.Nrows(), 1);
	ASSERT_EQ(mean.Ncols(), n_voxels);
	for (int i = 0; i < n_voxels; i++)
	{
		// GTEST has difficulty with comparing floats to 0
		ASSERT_TRUE(FloatEq(1, mean(1, i + 1) + 1));
	}

	mean = rundata.GetVoxelData("mean_c2");
	ASSERT_EQ(mean.Nrows(), 1);
	ASSERT_EQ(mean.Ncols(), n_voxels);
	for (int i = 0; i < n_voxels; i++)
	{
		ASSERT_TRUE(FloatEq(VAL * 1.5, mean(1, i + 1)));
	}

	mean = rundata.GetVoxelData("mean_c3");
	ASSERT_EQ(mean.Nrows(), 1);
	ASSERT_EQ(mean.Ncols(), n_voxels);
	for (int i = 0; i < n_voxels; i++)
	{
		ASSERT_TRUE(FloatEq(-VAL * 2, mean(1, i + 1)));
	}
}

// Test saving model fit data
TEST_P(InferenceMethodTest, SaveModelFit)
{
	int NTIMES = 10;
	int VSIZE = 5;
	float VAL = 2;
	int DEGREE = 3;
	int n_voxels = VSIZE * VSIZE * VSIZE;

	// Create coordinates and data matrices
	// Data fitted to a cubic function
	NEWMAT::Matrix voxelCoords, data;
	data.ReSize(NTIMES, n_voxels);
	voxelCoords.ReSize(3, n_voxels);
	int v = 1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				for (int n = 0; n < NTIMES; n++)
				{
					// function is VAL + 1.5VAL x n^2 - 2*VAL*n^3
					data(n + 1, v) = VAL + (1.5 * VAL) * (n + 1) * (n + 1) - 2 * VAL * (n + 1) * (n + 1) * (n + 1);
				}
				v++;
			}
		}
	}

	// Do just 1 iteration

	FabberIoMemory io;
	FabberRunData rundata(&io);
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data", data);
	rundata.Set("noise", "white");
	rundata.Set("model", "poly");
	rundata.Set("degree", stringify(DEGREE));
	rundata.Set("max-iterations", "10");
	rundata.SetBool("save-model-fit");
	rundata.Set("method", GetParam());
	rundata.Run();

	NEWMAT::Matrix mean = rundata.GetVoxelData("modelfit");
	ASSERT_EQ(mean.Nrows(), NTIMES);
	ASSERT_EQ(mean.Ncols(), n_voxels);
}

INSTANTIATE_TEST_CASE_P(MethodTests, InferenceMethodTest, ::testing::Values("vb", "nlls", "spatialvb"));

} // namespace

