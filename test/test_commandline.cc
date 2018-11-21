#ifndef NO_NEWIMAGE
// Tests fabber when run from the command line

#include "easylog.h"
#include "newimage/newimageall.h"
#include "rundata.h"
#include "setup.h"
#include "gtest/gtest.h"

#define ALLOWED_DELTA 0.001

namespace
{
// The fixture for testing class Foo.
class ClTestTest : public ::testing::TestWithParam<string>
{
protected:
    ClTestTest()
    {
    }

    virtual ~ClTestTest()
    {
    }

    virtual void SetUp()
    {
    }

    virtual void TearDown()
    {
        remove("testfabber.tmp");
        remove("out.tmp_latest");
        remove("uname.txt");
        int ret = system("rm -rf out.tmp"); // Not portable!
        if (ret != 0)
            cerr << "WARNING: failed to remove temp directory" << endl;
    }

    int runFabber(string cline)
    {
        string cmd = string(FABBER_BUILD_DIR) + "/fabber " + cline + ">testfabber.tmp";
        return system(cmd.c_str());
    }

    string getStdout()
    {
        ifstream in("testfabber.tmp");
        string str;

        str.assign((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        return str;
    }

    string getLogfile(string dir)
    {
        string file = dir + "/logfile";
        ifstream in(file.c_str());
        string str;

        str.assign((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        return str;
    }

    bool contains(string s, string txt)
    {
        return s.find(txt) != s.npos;
    }

    void compareNifti(string f1, string f2)
    {
        NEWIMAGE::volume<float> d1;
        read_volume(d1, f1);
        NEWIMAGE::volume<float> d2;
        read_volume(d2, f2);

        ASSERT_EQ(d1.xsize(), d2.xsize());
        ASSERT_EQ(d1.ysize(), d2.ysize());
        ASSERT_EQ(d1.zsize(), d2.zsize());

        //		Matrix m1 = d1.matrix();
        //		Matrix m2 = d2.matrix();

        for (int x = 1; x <= d1.xsize(); x++)
        {
            for (int y = 1; y <= d1.ysize(); y++)
            {
                for (int z = 1; z <= d1.zsize(); z++)
                {
                    ASSERT_NEAR(d1(x, y, z), d2(x, y, z), ALLOWED_DELTA);
                }
            }
        }
    }
};

// Tests help option
TEST_F(ClTestTest, Help)
{
    string args = "--help";

    ASSERT_EQ(0, runFabber(args));
    string out = getStdout();
    ASSERT_TRUE(contains(out, "Usage"));
}

// Test real data with a linear model and a VEST format basis
TEST_P(ClTestTest, LinearModelVest)
{
    string args = "--output=out.tmp --model=linear --basis=" + string(FABBER_SRC_DIR)
        + "/test/test_linear_design.mat ";
    args += " --mask=" + string(FABBER_SRC_DIR) + "/test/test_mask_small.nii.gz --data="
        + string(FABBER_SRC_DIR) + "/test/test_data.nii.gz --noise=white";
    args += " --method=" + GetParam() + " ";

    ASSERT_EQ(0, runFabber(args));
    string out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "model=linear"));
    ASSERT_TRUE(contains(out, "method=" + GetParam()));
    ASSERT_TRUE(contains(out, "test_mask_small.nii.gz"));
    ASSERT_TRUE(contains(out, "test_data.nii.gz"));

    compareNifti("out.tmp/mean_Parameter_1.nii.gz",
        string(FABBER_SRC_DIR) + "/test/outdata_linear_" + GetParam() + "/mean_Parameter_1.nii.gz");
    compareNifti("out.tmp/mean_Parameter_2.nii.gz",
        string(FABBER_SRC_DIR) + "/test/outdata_linear_" + GetParam() + "/mean_Parameter_2.nii.gz");
    compareNifti("out.tmp/mean_Parameter_3.nii.gz",
        string(FABBER_SRC_DIR) + "/test/outdata_linear_" + GetParam() + "/mean_Parameter_3.nii.gz");
    compareNifti("out.tmp/mean_Parameter_4.nii.gz",
        string(FABBER_SRC_DIR) + "/test/outdata_linear_" + GetParam() + "/mean_Parameter_4.nii.gz");
    compareNifti("out.tmp/zstat_Parameter_1.nii.gz", string(FABBER_SRC_DIR)
            + "/test/outdata_linear_" + GetParam() + "/zstat_Parameter_1.nii.gz");
    compareNifti("out.tmp/zstat_Parameter_2.nii.gz", string(FABBER_SRC_DIR)
            + "/test/outdata_linear_" + GetParam() + "/zstat_Parameter_2.nii.gz");
    compareNifti("out.tmp/zstat_Parameter_3.nii.gz", string(FABBER_SRC_DIR)
            + "/test/outdata_linear_" + GetParam() + "/zstat_Parameter_3.nii.gz");
    compareNifti("out.tmp/zstat_Parameter_4.nii.gz", string(FABBER_SRC_DIR)
            + "/test/outdata_linear_" + GetParam() + "/zstat_Parameter_4.nii.gz");
}

// Test real data with a linear model and an ASCII format basis
TEST_P(ClTestTest, LinearModelAscii)
{
    string args = "--output=out.tmp --model=linear --basis=" + string(FABBER_SRC_DIR)
        + "/test/test_linear_design_ascii.mat ";
    args += " --mask=" + string(FABBER_SRC_DIR) + "/test/test_mask_small.nii.gz --data="
        + string(FABBER_SRC_DIR) + "/test/test_data.nii.gz --noise=white";
    args += " --method=" + GetParam() + " ";

    ASSERT_EQ(0, runFabber(args));
    string out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "model=linear"));
    ASSERT_TRUE(contains(out, "method=" + GetParam()));
    ASSERT_TRUE(contains(out, "test_mask_small.nii.gz"));
    ASSERT_TRUE(contains(out, "test_data.nii.gz"));

    compareNifti("out.tmp/mean_Parameter_1.nii.gz",
        string(FABBER_SRC_DIR) + "/test/outdata_linear_" + GetParam() + "/mean_Parameter_1.nii.gz");
    compareNifti("out.tmp/mean_Parameter_2.nii.gz",
        string(FABBER_SRC_DIR) + "/test/outdata_linear_" + GetParam() + "/mean_Parameter_2.nii.gz");
    compareNifti("out.tmp/mean_Parameter_3.nii.gz",
        string(FABBER_SRC_DIR) + "/test/outdata_linear_" + GetParam() + "/mean_Parameter_3.nii.gz");
    compareNifti("out.tmp/mean_Parameter_4.nii.gz",
        string(FABBER_SRC_DIR) + "/test/outdata_linear_" + GetParam() + "/mean_Parameter_4.nii.gz");
    compareNifti("out.tmp/zstat_Parameter_1.nii.gz", string(FABBER_SRC_DIR)
            + "/test/outdata_linear_" + GetParam() + "/zstat_Parameter_1.nii.gz");
    compareNifti("out.tmp/zstat_Parameter_2.nii.gz", string(FABBER_SRC_DIR)
            + "/test/outdata_linear_" + GetParam() + "/zstat_Parameter_2.nii.gz");
    compareNifti("out.tmp/zstat_Parameter_3.nii.gz", string(FABBER_SRC_DIR)
            + "/test/outdata_linear_" + GetParam() + "/zstat_Parameter_3.nii.gz");
    compareNifti("out.tmp/zstat_Parameter_4.nii.gz", string(FABBER_SRC_DIR)
            + "/test/outdata_linear_" + GetParam() + "/zstat_Parameter_4.nii.gz");
}
// Test real data with a simple polynomial model
TEST_F(ClTestTest, PolyModel)
{
    string args = "--model=poly --output=out.tmp  --degree=2 --method=vb --noise=white ";
    args += " --mask=" + string(FABBER_SRC_DIR) + "/test/test_mask_small.nii.gz --data="
        + string(FABBER_SRC_DIR) + "/test/test_data.nii.gz";

    ASSERT_EQ(0, runFabber(args));
    string out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "model=poly"));
    ASSERT_TRUE(contains(out, "method=vb"));
    ASSERT_TRUE(contains(out, "noise=white"));
    ASSERT_TRUE(contains(out, "test_mask_small.nii.gz"));
    ASSERT_TRUE(contains(out, "test_data.nii.gz"));

    compareNifti(
        "out.tmp/mean_c0.nii.gz", string(FABBER_SRC_DIR) + "/test/outdata_poly/mean_c0.nii.gz");
    compareNifti(
        "out.tmp/mean_c1.nii.gz", string(FABBER_SRC_DIR) + "/test/outdata_poly/mean_c1.nii.gz");
    compareNifti(
        "out.tmp/mean_c2.nii.gz", string(FABBER_SRC_DIR) + "/test/outdata_poly/mean_c2.nii.gz");
    compareNifti(
        "out.tmp/std_c0.nii.gz", string(FABBER_SRC_DIR) + "/test/outdata_poly/std_c0.nii.gz");
    compareNifti(
        "out.tmp/std_c1.nii.gz", string(FABBER_SRC_DIR) + "/test/outdata_poly/std_c1.nii.gz");
    compareNifti(
        "out.tmp/std_c2.nii.gz", string(FABBER_SRC_DIR) + "/test/outdata_poly/std_c2.nii.gz");
}

// Test real data with an empty mask. May occur when running in parallel if some chunks have not ROI
TEST_F(ClTestTest, EmptyMask)
{
    string args = "--model=poly --output=out.tmp  --degree=2 --method=vb --noise=white ";
    args += " --mask=" + string(FABBER_SRC_DIR) + "/test/test_mask_empty.nii.gz --data="
        + string(FABBER_SRC_DIR) + "/test/test_data.nii.gz";

    ASSERT_EQ(0, runFabber(args));
    string out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "model=poly"));
    ASSERT_TRUE(contains(out, "method=vb"));
    ASSERT_TRUE(contains(out, "noise=white"));
    ASSERT_TRUE(contains(out, "test_mask_empty.nii.gz"));
    ASSERT_TRUE(contains(out, "test_data.nii.gz"));
}

// Test if output data retains properties from input
TEST_F(ClTestTest, OutputCopiesPropsMask)
{
    string args = "--model=poly --output=out.tmp  --degree=2 --method=vb --noise=white ";
    args += " --mask=" + string(FABBER_SRC_DIR) + "/test/test_mask_small.nii.gz --data="
        + string(FABBER_SRC_DIR) + "/test/test_data.nii.gz";

    ASSERT_EQ(0, runFabber(args));

    NEWIMAGE::volume<float> d1;
    read_volume(d1, string(FABBER_SRC_DIR) + "/test/test_data.nii.gz");
    NEWIMAGE::volume<float> d2;
    read_volume(d2, "out.tmp/mean_c0.nii.gz");
    ASSERT_EQ(d1.xdim(), d2.xdim());
}

// Test if output data retains properties from input
TEST_F(ClTestTest, OutputCopiesPropsNoMask)
{
    string args = "--model=poly --output=out.tmp  --degree=2 --method=vb --noise=white ";
    args += " --data=" + string(FABBER_SRC_DIR) + "/test/test_data_small.nii.gz";

    ASSERT_EQ(0, runFabber(args));

    NEWIMAGE::volume<float> d1;
    read_volume(d1, string(FABBER_SRC_DIR) + "/test/test_data_small.nii.gz");
    NEWIMAGE::volume<float> d2;
    read_volume(d2, "out.tmp/mean_c0.nii.gz");
    ASSERT_EQ(d1.xdim(), d2.xdim());
}

// Test fabber will run without a mask
TEST_F(ClTestTest, PolyModelNoMask)
{
    string args = "--model=poly --output=out.tmp  --degree=2 --method=vb --noise=white ";
    args += " --data=" + string(FABBER_SRC_DIR) + "/test/test_data_small.nii.gz";

    ASSERT_EQ(0, runFabber(args));
    string out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "model=poly"));
    ASSERT_TRUE(contains(out, "method=vb"));
    ASSERT_TRUE(contains(out, "noise=white"));
    ASSERT_TRUE(contains(out, "test_data_small.nii.gz"));
}

// Test fabber will append a + to the output dir if it already exists
TEST_F(ClTestTest, NoOverwrite)
{
    string args = "--model=poly --output=out.tmp  --degree=2 --method=vb --noise=white ";
    args += " --data=" + string(FABBER_SRC_DIR) + "/test/test_data_small.nii.gz";

    ASSERT_EQ(0, runFabber(args));
    string out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "model=poly"));
    ASSERT_TRUE(contains(out, "degree=2"));
    ASSERT_TRUE(contains(out, "method=vb"));
    ASSERT_TRUE(contains(out, "noise=white"));
    ASSERT_TRUE(contains(out, "test_data_small.nii.gz"));

    args = "--model=poly --output=out.tmp  --degree=1 --method=vb --noise=white ";
    args += " --data=" + string(FABBER_SRC_DIR) + "/test/test_data_small.nii.gz";

    ASSERT_EQ(0, runFabber(args));
    out = getLogfile("out.tmp+");
    ASSERT_TRUE(contains(out, "model=poly"));
    ASSERT_TRUE(contains(out, "degree=1"));
    ASSERT_TRUE(contains(out, "method=vb"));
    ASSERT_TRUE(contains(out, "noise=white"));
    ASSERT_TRUE(contains(out, "test_data_small.nii.gz"));

    int ret = system("rm -rf out.tmp+"); // Not portable!
    if (ret != 0)
        cerr << "WARNING: Failed to remove temp directory out.tmp+" << endl;
}

// Test fabber will overwrite the output dir if it already exists with --overwrite
TEST_F(ClTestTest, Overwrite)
{
    string args = "--model=poly --output=out.tmp  --degree=2 --method=vb --noise=white ";
    args += " --data=" + string(FABBER_SRC_DIR) + "/test/test_data_small.nii.gz";

    ASSERT_EQ(0, runFabber(args));
    string out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "model=poly"));
    ASSERT_TRUE(contains(out, "degree=2"));
    ASSERT_TRUE(contains(out, "method=vb"));
    ASSERT_TRUE(contains(out, "noise=white"));
    ASSERT_TRUE(contains(out, "test_data_small.nii.gz"));

    args = "--model=poly --output=out.tmp  --degree=1 --method=vb --noise=white ";
    args += " --data=" + string(FABBER_SRC_DIR) + "/test/test_data_small.nii.gz --overwrite";

    ASSERT_EQ(0, runFabber(args));
    out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "model=poly"));
    ASSERT_TRUE(contains(out, "degree=1"));
    ASSERT_TRUE(contains(out, "method=vb"));
    ASSERT_TRUE(contains(out, "noise=white"));
    ASSERT_TRUE(contains(out, "test_data_small.nii.gz"));
}

TEST_F(ClTestTest, ListModels)
{
    string args = "--listmodels";

    ASSERT_EQ(0, runFabber(args));
    string out = getStdout();
    ASSERT_TRUE(contains(out, "poly"));
    ASSERT_TRUE(contains(out, "linear"));
}

TEST_F(ClTestTest, ListMethods)
{
    string args = "--listmethods";

    ASSERT_EQ(0, runFabber(args));
    string out = getStdout();
    ASSERT_TRUE(contains(out, "nlls"));
    ASSERT_TRUE(contains(out, "vb"));
    ASSERT_TRUE(contains(out, "spatialvb"));
}

TEST_F(ClTestTest, UnusedParams)
{
    string args = "--model=poly --output=out.tmp  --degree=2 --method=vb --noise=white --overwrite";
    args += " --data=" + string(FABBER_SRC_DIR) + "/test/test_data_small.nii.gz";

    ASSERT_EQ(0, runFabber(args));
    string out = getLogfile("out.tmp");
    ASSERT_FALSE(contains(out, "WARNING"));

    args += " --squaffle";
    ASSERT_EQ(0, runFabber(args));
    out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "WARNING"));
    ASSERT_TRUE(contains(out, "Unused option"));
}

TEST_F(ClTestTest, OptFileOldStyle)
{
    std::ofstream outfile;
    outfile.open("optfile.tmp");
    outfile << "--model=poly" << endl
            << "--output=out.tmp" << endl
            << "--degree=2" << endl
            << "--method=vb" << endl
            << "--noise=white" << endl
            << "--overwrite" << endl;
    outfile.close();

    string args
        = "--optfile=optfile.tmp --data=" + string(FABBER_SRC_DIR) + "/test/test_data_small.nii.gz";

    ASSERT_EQ(0, runFabber(args));
    string out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "model=poly"));
    ASSERT_TRUE(contains(out, "degree=2"));
    ASSERT_TRUE(contains(out, "method=vb"));
    ASSERT_TRUE(contains(out, "noise=white"));
    ASSERT_TRUE(contains(out, "overwrite"));
    ASSERT_TRUE(contains(out, "test_data_small.nii.gz"));
}

TEST_F(ClTestTest, OptFileOldStyleOldName)
{
    std::ofstream outfile;
    outfile.open("optfile.tmp");
    outfile << "--model=poly" << endl
            << "--output=out.tmp" << endl
            << "--degree=2" << endl
            << "--method=vb" << endl
            << "--noise=white" << endl
            << "--overwrite" << endl;
    outfile.close();

    string args
        = "-@ optfile.tmp --data=" + string(FABBER_SRC_DIR) + "/test/test_data_small.nii.gz";

    ASSERT_EQ(0, runFabber(args));
    string out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "model=poly"));
    ASSERT_TRUE(contains(out, "degree=2"));
    ASSERT_TRUE(contains(out, "method=vb"));
    ASSERT_TRUE(contains(out, "noise=white"));
    ASSERT_TRUE(contains(out, "overwrite"));
    ASSERT_TRUE(contains(out, "test_data_small.nii.gz"));
}

TEST_F(ClTestTest, OptFileNewStyle)
{
    std::ofstream outfile;
    outfile.open("optfile.tmp");
    outfile << "model=poly" << endl
            << "output=out.tmp" << endl
            << "degree=2" << endl
            << "method=vb" << endl
            << "noise=white" << endl
            << "overwrite" << endl;
    outfile.close();

    string args
        = "-f optfile.tmp --data=" + string(FABBER_SRC_DIR) + "/test/test_data_small.nii.gz";

    ASSERT_EQ(0, runFabber(args));
    string out = getLogfile("out.tmp");
    ASSERT_TRUE(contains(out, "model=poly"));
    ASSERT_TRUE(contains(out, "degree=2"));
    ASSERT_TRUE(contains(out, "method=vb"));
    ASSERT_TRUE(contains(out, "noise=white"));
    ASSERT_TRUE(contains(out, "overwrite"));
    ASSERT_TRUE(contains(out, "test_data_small.nii.gz"));
}

INSTANTIATE_TEST_CASE_P(ClMethodTests, ClTestTest, ::testing::Values("vb", "nlls", "spatialvb"));
}
#endif
