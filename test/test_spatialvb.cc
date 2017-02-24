//
// Tests specifically for the Spatial VB method

#include "gtest/gtest.h"

#include "inference.h"
#include "inference_spatialvb.h"
#include "dataset.h"
#include "setup.h"

namespace
{

// Hack to allow us to get at private/protected members
class PublicVersion: public SpatialVariationalBayes
{
public:
	using SpatialVariationalBayes::CalcNeighbours;
	using SpatialVariationalBayes::m_neighbours;
	using SpatialVariationalBayes::m_neighbours2;
};

// The fixture for testing class Foo.
class SpatialVbTest: public ::testing::Test
{
protected:
	SpatialVbTest()
	{
		FabberSetup::SetupDefaults();
	}

	virtual ~SpatialVbTest()
	{
		FabberSetup::Destroy();
	}

	virtual void SetUp()
	{
		svb = reinterpret_cast<PublicVersion*> (static_cast<SpatialVariationalBayes*> (InferenceTechnique::NewFromName(
				"spatialvb")));
		model = FwdModel::NewFromName("poly");
		io.ClearVoxelData();
		rundata = new FabberRunData(&io);
		rundata->SetLogger(&log);
	}

	virtual void TearDown()
	{
		delete svb;
		delete model;
		delete rundata;
	}

	void Run()
	{
		io.Initialize(*rundata);

		std::auto_ptr<FwdModel> fwd_model(FwdModel::NewFromName(rundata->GetString("model")));
		fwd_model->Initialize(*rundata);

		svb->Initialize(fwd_model.get(), *rundata);
		svb->DoCalculations(*rundata);
		svb->SaveResults(*rundata);
	}

	void Initialize()
	{
		io.SetVoxelCoords(voxelCoords);
		rundata->Set("degree", "0");
		rundata->Set("noise", "white");
		svb->Initialize(model, *rundata);
	}

	FabberIoMemory io;
	NEWMAT::Matrix voxelCoords;
	FabberRunData *rundata;
	FwdModel *model;
	PublicVersion *svb;
	EasyLog log;
};

// Tests the CalcNeighbours method for a single voxel
TEST_F(SpatialVbTest, CalcNeighboursOneVoxel)
{
	// Create coordinates and data matrices
	voxelCoords.ReSize(3, 1);
	voxelCoords << 1 << 1 << 1;
	Initialize();

	svb->CalcNeighbours(voxelCoords);
	ASSERT_EQ(svb->m_neighbours.size(), 1);
	ASSERT_EQ(svb->m_neighbours[0].size(), 0);
}

// Tests the CalcNeighbours method for a single voxel at zero
TEST_F(SpatialVbTest, CalcNeighboursOneVoxelZero)
{
	// Create coordinates and data matrices
	voxelCoords.ReSize(3, 1);
	voxelCoords << 0 << 0 << 0;
	Initialize();

	svb->CalcNeighbours(voxelCoords);
	ASSERT_EQ(svb->m_neighbours.size(), 1);
	ASSERT_EQ(svb->m_neighbours[0].size(), 0);
}

// Tests the CalcNeighbours method for multi voxels in X direction
TEST_F(SpatialVbTest, CalcNeighboursMultiVoxelsX)
{
	int NVOXELS = 5;

	// Create coordinates and data matrices
	voxelCoords.ReSize(3, NVOXELS);
	for (int v=1; v<=NVOXELS; v++)
	{
		voxelCoords(1, v) = v;
		voxelCoords(2, v) = 1;
		voxelCoords(3, v) = 1;
	}
	Initialize();

	svb->CalcNeighbours(voxelCoords);
	ASSERT_EQ(svb->m_neighbours.size(), NVOXELS);

	for (int v=1; v<=NVOXELS; v++)
	{
		int expected=2;
		if (v == 1 || v == NVOXELS) expected=1;
		ASSERT_EQ(svb->m_neighbours[v-1].size(), expected);
	}
}

// Tests the CalcNeighbours method for multi voxels in Y direction
TEST_F(SpatialVbTest, CalcNeighboursMultiVoxelsY)
{
	int NVOXELS = 5;

	// Create coordinates and data matrices
	voxelCoords.ReSize(3, NVOXELS);
	for (int v=1; v<=NVOXELS; v++)
	{
		voxelCoords(1, v) = 1;
		voxelCoords(2, v) = v;
		voxelCoords(3, v) = 1;
	}
	Initialize();

	svb->CalcNeighbours(voxelCoords);
	ASSERT_EQ(svb->m_neighbours.size(), NVOXELS);

	for (int v=1; v<=NVOXELS; v++)
	{
		int expected=2;
		if (v == 1 || v == NVOXELS) expected=1;
		ASSERT_EQ(svb->m_neighbours[v-1].size(), expected);
	}
}

// Tests the CalcNeighbours method for multi voxels in Z direction
TEST_F(SpatialVbTest, CalcNeighboursMultiVoxelsZ)
{
	int NVOXELS = 5;

	// Create coordinates and data matrices
	voxelCoords.ReSize(3, NVOXELS);
	for (int v=1; v<=NVOXELS; v++)
	{
		voxelCoords(1, v) = 1;
		voxelCoords(2, v) = 1;
		voxelCoords(3, v) = v;
	}
	Initialize();

	svb->CalcNeighbours(voxelCoords);
	ASSERT_EQ(svb->m_neighbours.size(), NVOXELS);

	for (int v=1; v<=NVOXELS; v++)
	{
		int expected=2;
		if (v == 1 || v == NVOXELS) expected=1;
		ASSERT_EQ(svb->m_neighbours[v-1].size(), expected);
	}
}

// Tests the CalcNeighbours method for a cubic volume
// with coords starting at zero
TEST_F(SpatialVbTest, CalcNeighboursMultiVoxels3dZeros)
{
	int VSIZE = 5;
	int NVOXELS = VSIZE*VSIZE*VSIZE;

	voxelCoords.ReSize(3, NVOXELS);
	int v=1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				v++;
			}
		}
	}

	Initialize();
	ASSERT_NO_THROW(svb->CalcNeighbours(voxelCoords));
	ASSERT_EQ(svb->m_neighbours.size(), NVOXELS);

	v=1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				vector<int> expected;
				if (x != 0) expected.push_back(v-1);
				if (x != VSIZE-1) expected.push_back(v+1);
				if (y != 0) expected.push_back(v-VSIZE);
				if (y != VSIZE-1) expected.push_back(v+VSIZE);
				if (z != 0) expected.push_back(v-VSIZE*VSIZE);
				if (z != VSIZE-1) expected.push_back(v+VSIZE*VSIZE);
				//                cout << "voxel " << v << " expected=" << expected.size() << endl;
				ASSERT_EQ(svb->m_neighbours[v-1].size(), expected.size());
				vector<int> nb = svb->m_neighbours[v-1];
				for (vector<int>::iterator iter=nb.begin(); iter!=nb.end(); iter++)
				{
					vector<int>::iterator found = std::find(expected.begin(), expected.end(), *iter);
					ASSERT_TRUE(found != expected.end());
				}
				v++;
			}
		}
	}
}

// Tests the CalcNeighbours method with reduced spatial
// dimensions (1d)
TEST_F(SpatialVbTest, CalcNeighboursCubic1D)
{
	int VSIZE = 5;
	int NVOXELS = VSIZE*VSIZE*VSIZE;

	voxelCoords.ReSize(3, NVOXELS);
	int v=1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				v++;
			}
		}
	}

	rundata->Set("spatial-dims", "1");
	Initialize();
	ASSERT_NO_THROW(svb->CalcNeighbours(voxelCoords));
	ASSERT_EQ(svb->m_neighbours.size(), NVOXELS);

	v=1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				vector<int> expected;
				if (x != 0) expected.push_back(v-1);
				if (x != VSIZE-1) expected.push_back(v+1);
				//                cout << "voxel " << v << " expected=" << expected.size() << endl;
				ASSERT_EQ(svb->m_neighbours[v-1].size(), expected.size());
				vector<int> nb = svb->m_neighbours[v-1];
				for (vector<int>::iterator iter=nb.begin(); iter!=nb.end(); iter++)
				{
					vector<int>::iterator found = std::find(expected.begin(), expected.end(), *iter);
					ASSERT_TRUE(found != expected.end());
				}
				v++;
			}
		}
	}
}

// Tests the CalcNeighbours method with reduced spatial
// dimensions (2d)
TEST_F(SpatialVbTest, CalcNeighboursCubic2D)
{
	int VSIZE = 5;
	int NVOXELS = VSIZE*VSIZE*VSIZE;

	voxelCoords.ReSize(3, NVOXELS);
	int v=1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				v++;
			}
		}
	}

	rundata->Set("spatial-dims", "2");
	Initialize();
	ASSERT_NO_THROW(svb->CalcNeighbours(voxelCoords));
	ASSERT_EQ(svb->m_neighbours.size(), NVOXELS);

	v=1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				vector<int> expected;
				if (x != 0) expected.push_back(v-1);
				if (x != VSIZE-1) expected.push_back(v+1);
				if (y != 0) expected.push_back(v-VSIZE);
				if (y != VSIZE-1) expected.push_back(v+VSIZE);
				//                cout << "voxel " << v << " expected=" << expected.size() << endl;
				ASSERT_EQ(svb->m_neighbours[v-1].size(), expected.size());
				vector<int> nb = svb->m_neighbours[v-1];
				for (vector<int>::iterator iter=nb.begin(); iter!=nb.end(); iter++)
				{
					vector<int>::iterator found = std::find(expected.begin(), expected.end(), *iter);
					ASSERT_TRUE(found != expected.end());
				}
				v++;
			}
		}
	}
}

// Tests the second-nearest-neighbour calculation
// for cubic cell
TEST_F(SpatialVbTest, CalcNeighbours2Cubic)
{
	int VSIZE = 5;
	int NVOXELS = VSIZE*VSIZE*VSIZE;

	voxelCoords.ReSize(3, NVOXELS);
	int v=1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x;
				voxelCoords(2, v) = y;
				voxelCoords(3, v) = z;
				v++;
			}
		}
	}

	Initialize();
	ASSERT_NO_THROW(svb->CalcNeighbours(voxelCoords));
	ASSERT_EQ(svb->m_neighbours.size(), NVOXELS);

	v=1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				vector<int> expected;
				if (x >= 2) expected.push_back(v-2);
				if (x <= VSIZE-3) expected.push_back(v+2);
				if (y >= 2) expected.push_back(v-(2*VSIZE));
				if (y <= VSIZE-3) expected.push_back(v+(2*VSIZE));
				if (z >=2) expected.push_back(v-(2*VSIZE*VSIZE));
				if (z <= VSIZE-3) expected.push_back(v+(2*VSIZE*VSIZE));

				// Neighbours on a diagonal. Note that all of these can be
				// reached in two ways and current code keeps these duplicates,
				// so we expect them here too.
				for (int i=0; i<2; i++)
				{
					if ((x>=1) && (y >= 1)) expected.push_back(v-1-VSIZE);
					if ((x>=1) && (y <= VSIZE-2)) expected.push_back(v-1+VSIZE);
					if ((x<= VSIZE-2) && (y >= 1)) expected.push_back(v+1-VSIZE);
					if ((x<= VSIZE-2) && (y <= VSIZE-2)) expected.push_back(v+1+VSIZE);

					if ((x>=1) && (z >= 1)) expected.push_back(v-1-(VSIZE*VSIZE));
					if ((x>=1) && (z <= VSIZE-2)) expected.push_back(v-1+(VSIZE*VSIZE));
					if ((x<= VSIZE-2) && (z >= 1)) expected.push_back(v+1-(VSIZE*VSIZE));
					if ((x<= VSIZE-2) && (z <= VSIZE-2)) expected.push_back(v+1+(VSIZE*VSIZE));

					if ((y>=1) && (z >= 1)) expected.push_back(v-VSIZE-(VSIZE*VSIZE));
					if ((y>=1) && (z <= VSIZE-2)) expected.push_back(v-VSIZE+(VSIZE*VSIZE));
					if ((y<= VSIZE-2) && (z >= 1)) expected.push_back(v+VSIZE-(VSIZE*VSIZE));
					if ((y<= VSIZE-2) && (z <= VSIZE-2)) expected.push_back(v+VSIZE+(VSIZE*VSIZE));
				}
				//                cout << "voxel " << v << " expected=" << expected.size() << endl;
				ASSERT_EQ(svb->m_neighbours2[v-1].size(), expected.size());
				vector<int> nb = svb->m_neighbours2[v-1];
				for (vector<int>::iterator iter=nb.begin(); iter!=nb.end(); iter++)
				{
					vector<int>::iterator found = std::find(expected.begin(), expected.end(), *iter);
					ASSERT_TRUE(found != expected.end());
				}
				v++;
			}
		}
	}
}

// Tests the CalcNeighbours method for a cubic volume
// with coords starting at 1
TEST_F(SpatialVbTest, CalcNeighboursMultiVoxels3dNoZeros)
{
	int VSIZE = 5;
	int NVOXELS = VSIZE*VSIZE*VSIZE;

	voxelCoords.ReSize(3, NVOXELS);
	int v=1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				voxelCoords(1, v) = x+1;
				voxelCoords(2, v) = y+1;
				voxelCoords(3, v) = z+1;
				v++;
			}
		}
	}

	Initialize();
	svb->CalcNeighbours(voxelCoords);
	ASSERT_EQ(svb->m_neighbours.size(), NVOXELS);

	v=1;
	for (int z = 0; z < VSIZE; z++)
	{
		for (int y = 0; y < VSIZE; y++)
		{
			for (int x = 0; x < VSIZE; x++)
			{
				vector<int> expected;
				if (x != 0) expected.push_back(v-1);
				if (x != VSIZE-1) expected.push_back(v+1);
				if (y != 0) expected.push_back(v-VSIZE);
				if (y != VSIZE-1) expected.push_back(v+VSIZE);
				if (z != 0) expected.push_back(v-VSIZE*VSIZE);
				if (z != VSIZE-1) expected.push_back(v+VSIZE*VSIZE);
				//                cout << "voxel " << v << " expected=" << expected.size() << endl;
				ASSERT_EQ(svb->m_neighbours[v-1].size(), expected.size());
				vector<int> nb = svb->m_neighbours[v-1];
				for (vector<int>::iterator iter=nb.begin(); iter!=nb.end(); iter++)
				{
					vector<int>::iterator found = std::find(expected.begin(), expected.end(), *iter);
					ASSERT_TRUE(found != expected.end());
				}
				v++;
			}
		}
	}
}

// Tests the CalcNeighbours method for an irregular volume
TEST_F(SpatialVbTest, CalcNeighboursMultiVoxels3dIrregular)
{
	int NVOXELS = 5;

	voxelCoords.ReSize(3, NVOXELS);
	voxelCoords << 1 << 2 << 1 << 2 << 1
	<< 1 << 1 << 2 << 2 << 1
	<< 1 << 1 << 1 << 1 << 2;

	Initialize();
	svb->CalcNeighbours(voxelCoords);
	ASSERT_EQ(svb->m_neighbours.size(), NVOXELS);

	vector<int> n;
	n = svb->m_neighbours[0];
	ASSERT_EQ(n.size(), 3);

	n = svb->m_neighbours[1];
	ASSERT_EQ(n.size(), 2);

	n = svb->m_neighbours[2];
	ASSERT_EQ(n.size(), 2);

	n = svb->m_neighbours[3];
	ASSERT_EQ(n.size(), 2);

	n = svb->m_neighbours[4];
	ASSERT_EQ(n.size(), 1);
}

// Tests the CalcNeighbours method for an irregular volume
// looking at second-neighbours
TEST_F(SpatialVbTest, CalcNeighbours2MultiVoxels3dIrregular)
{
	int NVOXELS = 5;

	voxelCoords.ReSize(3, NVOXELS);
	voxelCoords
	<< 1 << 2 << 1 << 2 << 1
	<< 1 << 1 << 2 << 2 << 1
	<< 1 << 1 << 1 << 1 << 2;

	Initialize();
	svb->CalcNeighbours(voxelCoords);
	ASSERT_EQ(svb->m_neighbours.size(), NVOXELS);

	vector<int> n;
	// 1 neighbour2, two ways to get there
	n = svb->m_neighbours2[0];
	ASSERT_EQ(n.size(), 2);

	// 2 neighbour2s, one has two ways to get there
	n = svb->m_neighbours2[1];
	ASSERT_EQ(n.size(), 3);

	// 2 neighbour2s, one has two ways to get there
	n = svb->m_neighbours2[2];
	ASSERT_EQ(n.size(), 3);

	// 1 neighbour2s, two ways to get there
	n = svb->m_neighbours2[3];
	ASSERT_EQ(n.size(), 2);

	// 2 neighbour2s, one way to get there
	n = svb->m_neighbours2[4];
	ASSERT_EQ(n.size(), 2);
}

// Test restarting VB run with the output-only option
TEST_F(SpatialVbTest, RestartOutputOnly)
{
	int NTIMES = 10; // needs to be even
	int VSIZE = 5;
	float VAL = 7.32;
	int REPEATS = 50;
	int DEGREE = 2;

	// Create coordinates and data matrices
	// Data fitted to a quadratic function
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
					data(n + 1, v) = VAL + (1.5 * VAL) * (n + 1) * (n + 1);
				}
				v++;
			}
		}
	}

	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data",data);
	rundata->Set("noise", "white");
	rundata->Set("model", "poly");
	rundata->Set("degree", stringify(DEGREE));
	rundata->Set("method", "spatialvb");
	rundata->Set("max-iterations", "50");
	rundata->SetBool("save-mean", false);
	rundata->SetBool("save-modelfit", false);
	rundata->SetBool("save-mvn", true);
	rundata->Run();

	NEWMAT::Matrix mvns = rundata->GetVoxelData("finalMVN");
	TearDown();
	SetUp();
	rundata->Set("noise", "white");
	rundata->Set("model", "poly");
	rundata->Set("degree", stringify(DEGREE));
	rundata->Set("method", "spatialvb");
	rundata->SetBool("save-mean", true);
	rundata->SetBool("save-model-fit", true);
	rundata->SetBool("save-mvn", false);
	rundata->SetBool("output-only", true);
	rundata->Set("continue-from-mvn", "mvns");
	io.SetVoxelCoords(voxelCoords);
	io.SetVoxelData("data",data);
	io.SetVoxelData("mvns", mvns);
	Run();

	NEWMAT::Matrix mean = rundata->GetVoxelData("mean_c0");
	ASSERT_EQ(mean.Nrows(), 1);
	ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
	for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
	{
		ASSERT_NEAR(VAL, mean(1, i + 1), 0.0001);
	}
	mean = rundata->GetVoxelData("mean_c1");
	ASSERT_EQ(mean.Nrows(), 1);
	ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
	for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
	{
		// GTEST has difficulty with comparing floats to 0
		ASSERT_NEAR(1, mean(1, i + 1) + 1, 0.0001);
	}
	mean = rundata->GetVoxelData("mean_c2");
	ASSERT_EQ(mean.Nrows(), 1);
	ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
	for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
	{
		ASSERT_NEAR(VAL * 1.5, mean(1, i + 1), 0.0001);
	}

	NEWMAT::Matrix fit = rundata->GetVoxelData("modelfit");
	ASSERT_EQ(fit.Nrows(), NTIMES);
	ASSERT_EQ(fit.Ncols(), VSIZE * VSIZE * VSIZE);
}

} // namespace

