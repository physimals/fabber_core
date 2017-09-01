//
// Tests of the run data class

#include "gtest/gtest.h"

#include "easylog.h"
#include "rundata.h"
#include "setup.h"

#include <fstream>

namespace
{
// The fixture for testing class Foo.
class RunDataTest : public ::testing::Test
{
protected:
    RunDataTest()
    {
    }

    virtual ~RunDataTest()
    {
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

// Tests concatenated data sets
TEST_F(RunDataTest, ConcatenatedData)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data1, data2, data3;
    data1.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    data2.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    data3.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
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
                    data1(n + 1, v) = VAL;
                    data2(n + 1, v) = VAL * 2;
                    data3(n + 1, v) = VAL * 3;
                }
                v++;
            }
        }
    }

    EasyLog log;
    FabberRunData rundata;
    rundata.SetLogger(&log);
    rundata.SetVoxelCoords(voxelCoords);
    rundata.SetVoxelData("data1", data1);
    rundata.SetVoxelData("data2", data2);
    rundata.SetVoxelData("data3", data3);
    rundata.Set("data-order", "concatenate");
    NEWMAT::Matrix data = rundata.GetMainVoxelData();

    ASSERT_EQ(data.Nrows(), NTIMES * 3);
    ASSERT_EQ(data.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        for (int t = 0; t < NTIMES * 3; t++)
        {
            if (t < NTIMES)
                ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL);
            else if (t < NTIMES * 2)
                ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL * 2);
            else
                ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL * 3);
        }
    }
}

// Tests interleaved data sets
TEST_F(RunDataTest, InterleavedData)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data1, data2, data3;
    data1.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    data2.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    data3.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
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
                    data1(n + 1, v) = VAL;
                    data2(n + 1, v) = VAL * 2;
                    data3(n + 1, v) = VAL * 3;
                }
                v++;
            }
        }
    }

    EasyLog log;
    FabberRunData rundata;
    rundata.SetLogger(&log);
    rundata.SetVoxelCoords(voxelCoords);
    rundata.SetVoxelData("data1", data1);
    rundata.SetVoxelData("data2", data2);
    rundata.SetVoxelData("data3", data3);
    rundata.Set("data-order", "interleave");
    NEWMAT::Matrix data = rundata.GetMainVoxelData();

    ASSERT_EQ(data.Nrows(), NTIMES * 3);
    ASSERT_EQ(data.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        for (int t = 0; t < NTIMES * 3; t++)
        {
            if (t % 3 == 0)
                ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL);
            else if (t % 3 == 1)
                ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL * 2);
            else
                ASSERT_FLOAT_EQ(data(t + 1, i + 1), VAL * 3);
        }
    }
}
// Tests inconsistent multi-data sets
TEST_F(RunDataTest, MultiDataInconsistent)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data1, data2, data3;
    data1.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    data2.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    data3.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
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
                    data1(n + 1, v) = VAL;
                    data2(n + 1, v) = VAL * 2;
                    data3(n + 1, v) = VAL * 3;
                }
                v++;
            }
        }
    }

    FabberRunData rundata;
    rundata.SetVoxelCoords(voxelCoords);
    rundata.SetVoxelData("data1", data1);
    rundata.SetVoxelData("data2", data2);
    rundata.SetVoxelData("data3", data3);
    rundata.Set("data-order", "singlefile");
    ASSERT_THROW(NEWMAT::Matrix data = rundata.GetMainVoxelData(), InvalidOptionValue);
}

// Tests that we can read config from a .fab file
TEST_F(RunDataTest, OptionsFile)
{
    string FILENAME = "test_config";

    ofstream os;
    os.open(FILENAME.c_str(), ios::out);
    os << "noise=white" << endl
       << "model=poly" << endl
       << "method=vb" << endl
       << "bool-option" << endl
       << "#comment, ignored" << endl;
    os.close();

    FabberRunData rundata;
    rundata.ParseParamFile(FILENAME);
    ASSERT_EQ("white", rundata.GetString("noise"));
    ASSERT_EQ("poly", rundata.GetString("model"));
    ASSERT_EQ("vb", rundata.GetString("method"));
    ASSERT_EQ(true, rundata.GetBool("bool-option"));
}

// Tests embedded comments
TEST_F(RunDataTest, OptionsFileEmbeddedComment)
{
    string FILENAME = "test_config";

    ofstream os;
    os.open(FILENAME.c_str(), ios::out);
    os << "model=poly" << endl
       << "degree=0 # Keep things simple" << endl;
    os.close();

    FabberRunData rundata;
    rundata.ParseParamFile(FILENAME);
    ASSERT_EQ("poly", rundata.GetString("model"));
}

// Tests unsetting an option
TEST_F(RunDataTest, Unset)
{
    FabberRunData rundata;
    rundata.Set("wibble", "wobble");
    rundata.SetBool("bobble");

    ASSERT_EQ("wobble", rundata.GetStringDefault("wibble", "squabble"));
    rundata.Unset("wibble");
    ASSERT_EQ("squabble", rundata.GetStringDefault("wibble", "squabble"));
    ASSERT_EQ(true, rundata.GetBool("bobble"));
    rundata.Unset("bobble");
    ASSERT_EQ(false, rundata.GetBool("bobble"));
}

// Tests integer option
TEST_F(RunDataTest, IntConvert)
{
    FabberRunData rundata;
    rundata.Set("wibble", "7");

    ASSERT_EQ(7, rundata.GetInt("wibble"));
    ASSERT_EQ(7, rundata.GetInt("wibble", 7, 8));
    ASSERT_EQ(7, rundata.GetInt("wibble", 6, 7));
    ASSERT_EQ(7, rundata.GetInt("wibble", 7, 7));
}

// Tests integer list option
TEST_F(RunDataTest, IntList)
{
    FabberRunData rundata;
    rundata.Set("wibble1", "7");
    rundata.Set("wibble2", "8");

    vector<int> list = rundata.GetIntList("wibble");
    ASSERT_EQ(2, list.size());
    ASSERT_EQ(7, list[0]);
    ASSERT_EQ(8, list[1]);
}

// Tests double option
TEST_F(RunDataTest, DoubleConvert)
{
    FabberRunData rundata;
    rundata.Set("wibble", "7.5");

    ASSERT_EQ(7.5, rundata.GetDouble("wibble"));
    ASSERT_EQ(7.5, rundata.GetDouble("wibble", 7, 8));
}

// Tests double list option
TEST_F(RunDataTest, DoubleList)
{
    FabberRunData rundata;
    rundata.Set("wibble1", "7.5");
    rundata.Set("wibble2", "8.5");

    vector<double> list = rundata.GetDoubleList("wibble");
    ASSERT_EQ(2, list.size());
    ASSERT_EQ(7.5, list[0]);
    ASSERT_EQ(8.5, list[1]);
}

// Tests bad int option
TEST_F(RunDataTest, IntConvertFail)
{
    FabberRunData rundata;
    rundata.Set("wibble", "ABC");

    ASSERT_THROW(rundata.GetInt("wibble"), InvalidOptionValue);
}

// Tests bad int list option
TEST_F(RunDataTest, IntListConvertFail)
{
    FabberRunData rundata;
    rundata.Set("wibble1", "2");
    rundata.Set("wibble2", "ABC");

    ASSERT_THROW(rundata.GetIntList("wibble"), InvalidOptionValue);
}

// Tests bad double option
TEST_F(RunDataTest, DoubleConvertFail)
{
    FabberRunData rundata;
    rundata.Set("wibble", "ABC");

    ASSERT_THROW(rundata.GetDouble("wibble"), InvalidOptionValue);
}

// Tests bad double list option
TEST_F(RunDataTest, DoubleListConvertFail)
{
    FabberRunData rundata;
    rundata.Set("wibble1", "2.3");
    rundata.Set("wibble2", "ABC");

    ASSERT_THROW(rundata.GetDoubleList("wibble"), InvalidOptionValue);
}

// Tests int option that's too small
TEST_F(RunDataTest, IntConvertFailMin)
{
    FabberRunData rundata;
    rundata.Set("wibble", "7");

    ASSERT_THROW(rundata.GetInt("wibble", 10), InvalidOptionValue);
}

// Tests int option that's too big
TEST_F(RunDataTest, IntConvertFailMax)
{
    FabberRunData rundata;
    rundata.Set("wibble", "7");

    ASSERT_THROW(rundata.GetInt("wibble", 0, 3), InvalidOptionValue);
}

// Tests int list option that's too small
TEST_F(RunDataTest, IntListConvertFailMin)
{
    FabberRunData rundata;
    rundata.Set("wibble1", "14");
    rundata.Set("wibble2", "7");

    ASSERT_THROW(rundata.GetIntList("wibble", 10), InvalidOptionValue);
}

// Tests int list option that's too big
TEST_F(RunDataTest, IntListConvertFailMax)
{
    FabberRunData rundata;
    rundata.Set("wibble1", "1");
    rundata.Set("wibble2", "2");
    rundata.Set("wibble3", "7");

    ASSERT_THROW(rundata.GetIntList("wibble", 0, 3), InvalidOptionValue);
}

// Tests double option that's too small
TEST_F(RunDataTest, DoubleConvertFailMin)
{
    FabberRunData rundata;
    rundata.Set("wibble", "7.5");

    ASSERT_THROW(rundata.GetDouble("wibble", 10.1), InvalidOptionValue);
}

// Tests double option that's too big
TEST_F(RunDataTest, DoubleConvertFailMax)
{
    FabberRunData rundata;
    rundata.Set("wibble", "7.5");

    ASSERT_THROW(rundata.GetDouble("wibble", 0.2, 3.3), InvalidOptionValue);
}

// Tests double list option that's too small
TEST_F(RunDataTest, DoubleListConvertFailMin)
{
    FabberRunData rundata;
    rundata.Set("wibble1", "14.8");
    rundata.Set("wibble2", "7.3");

    ASSERT_THROW(rundata.GetDoubleList("wibble", 10), InvalidOptionValue);
}

// Tests double list option that's too big
TEST_F(RunDataTest, DoubleListConvertFailMax)
{
    FabberRunData rundata;
    rundata.Set("wibble1", "1.2");
    rundata.Set("wibble2", "2.3");
    rundata.Set("wibble3", "7.1");

    ASSERT_THROW(rundata.GetDoubleList("wibble", 0, 3), InvalidOptionValue);
}

// Tests odd case where data name could lead to circular reference
TEST_F(RunDataTest, CircularDataRef)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data1;
    data1.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
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
                    data1(n + 1, v) = VAL;
                }
                v++;
            }
        }
    }

    FabberRunData rundata;
    rundata.SetVoxelCoords(voxelCoords);
    rundata.SetVoxelData("data", data1);
    rundata.Set("data", "data");

    NEWMAT::Matrix data = rundata.GetMainVoxelData();

    ASSERT_EQ(VSIZE * VSIZE * VSIZE, data.Ncols());
    ASSERT_EQ(NTIMES, data.Nrows());
}

// Tests clearing voxel data
TEST_F(RunDataTest, ClearVoxelData)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data1, data2, data3;
    data1.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    data2.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    data3.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
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
                    data1(n + 1, v) = VAL;
                    data2(n + 1, v) = VAL * 2;
                    data3(n + 1, v) = VAL * 3;
                }
                v++;
            }
        }
    }

    FabberRunData rundata;
    rundata.SetVoxelCoords(voxelCoords);
    rundata.SetVoxelData("data1", data1);
    rundata.SetVoxelData("data2", data2);
    rundata.SetVoxelData("data3", data3);

    rundata.ClearVoxelData("data1");
    ASSERT_NO_THROW(rundata.GetVoxelCoords());
    ASSERT_THROW(rundata.GetVoxelData("data1"), DataNotFound);
    ASSERT_NO_THROW(rundata.GetVoxelData("data2"));
    ASSERT_NO_THROW(rundata.GetVoxelData("data3"));

    rundata.ClearVoxelData();
    ASSERT_THROW(rundata.GetVoxelCoords(), DataNotFound);
    ASSERT_THROW(rundata.GetVoxelData("data1"), DataNotFound);
    ASSERT_THROW(rundata.GetVoxelData("data2"), DataNotFound);
    ASSERT_THROW(rundata.GetVoxelData("data3"), DataNotFound);
}
}
