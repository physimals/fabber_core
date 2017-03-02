#include "gtest/gtest.h"

#include "inference_vb.h"
#include "rundata.h"

#include "newmat.h"

class PriorsTest : public ::testing::Test {
protected:
    PriorsTest()
    {
    }

    virtual ~PriorsTest()
    {
    }

    virtual void SetUp()
    {
    }

    virtual void TearDown()
    {
    }
};

TEST_F(PriorsTest, DefaultPriors)
{
    string PARAM_NAME = "abc";
    int PARAM_IDX = 7;

    FabberRunData rundata(0);

    PriorType prior(PARAM_IDX, PARAM_NAME, rundata);
    ASSERT_EQ(PARAM_NAME, prior.m_param_name);
    ASSERT_EQ(PARAM_IDX, prior.m_idx);
    ASSERT_EQ('-', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("", prior.m_filename);
    ASSERT_EQ(0, prior.m_image.Ncols());
}

TEST_F(PriorsTest, MultiplePriorsSpatialMethod)
{
    FabberRunData rundata(0);
    rundata.Set("param-spatial-priors", "ABC");

    PriorType prior(0, "a", rundata);
    ASSERT_EQ("a", prior.m_param_name);
    ASSERT_EQ(0, prior.m_idx);
    ASSERT_EQ('A', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("", prior.m_filename);
    ASSERT_EQ(0, prior.m_image.Ncols());

    prior = PriorType(1, "b", rundata);
    ASSERT_EQ("b", prior.m_param_name);
    ASSERT_EQ(1, prior.m_idx);
    ASSERT_EQ('B', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("", prior.m_filename);
    ASSERT_EQ(0, prior.m_image.Ncols());

    prior = PriorType(2, "c", rundata);
    ASSERT_EQ("c", prior.m_param_name);
    ASSERT_EQ(2, prior.m_idx);
    ASSERT_EQ('C', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("", prior.m_filename);
    ASSERT_EQ(0, prior.m_image.Ncols());
}

TEST_F(PriorsTest, SinglePlus)
{
    FabberRunData rundata(0);
    rundata.Set("param-spatial-priors", "+");

    PriorType prior(0, "a", rundata);
    ASSERT_EQ("a", prior.m_param_name);
    ASSERT_EQ(0, prior.m_idx);
    ASSERT_EQ('-', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("", prior.m_filename);
}

TEST_F(PriorsTest, MultiplePriorsSpatialMethodPlus)
{
    FabberRunData rundata(0);
    rundata.Set("param-spatial-priors", "AB+");

    PriorType prior(0, "a", rundata);
    ASSERT_EQ("a", prior.m_param_name);
    ASSERT_EQ(0, prior.m_idx);
    ASSERT_EQ('A', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("", prior.m_filename);
    ASSERT_EQ(0, prior.m_image.Ncols());

    prior = PriorType(1, "b", rundata);
    ASSERT_EQ("b", prior.m_param_name);
    ASSERT_EQ(1, prior.m_idx);
    ASSERT_EQ('B', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("", prior.m_filename);
    ASSERT_EQ(0, prior.m_image.Ncols());

    prior = PriorType(2, "c", rundata);
    ASSERT_EQ("c", prior.m_param_name);
    ASSERT_EQ(2, prior.m_idx);
    ASSERT_EQ('B', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("", prior.m_filename);
    ASSERT_EQ(0, prior.m_image.Ncols());

    prior = PriorType(3, "d", rundata);
    ASSERT_EQ("d", prior.m_param_name);
    ASSERT_EQ(3, prior.m_idx);
    ASSERT_EQ('B', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("", prior.m_filename);
}

TEST_F(PriorsTest, ImagePriorSpatialMethod)
{
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix iprior_data1, iprior_data2;
    iprior_data1.ReSize(1, VSIZE * VSIZE * VSIZE);
    iprior_data2.ReSize(1, VSIZE * VSIZE * VSIZE);
    int v = 1;
    for (int z = 0; z < VSIZE; z++) {
        for (int y = 0; y < VSIZE; y++) {
            for (int x = 0; x < VSIZE; x++) {
                iprior_data1(1, v) = VAL * 1.5;
                iprior_data2(1, v) = VAL * 2.5;
                v++;
            }
        }
    }

    FabberRunData rundata;
    rundata.SetVoxelData("PSP_byname1_image", iprior_data1);
    rundata.Set("param-spatial-priors", "I");

    PriorType prior(0, "a", rundata);
    ASSERT_EQ("a", prior.m_param_name);
    ASSERT_EQ(0, prior.m_idx);
    ASSERT_EQ('I', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("PSP_byname1_image", prior.m_filename);
    NEWMAT::RowVector image = prior.m_image;
    ASSERT_EQ(VSIZE * VSIZE * VSIZE, image.Ncols());
}

TEST_F(PriorsTest, ImagePriorSpatialMethodMultiple)
{
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix iprior_data1, iprior_data2;
    iprior_data1.ReSize(1, VSIZE * VSIZE * VSIZE);
    iprior_data2.ReSize(1, VSIZE * VSIZE * VSIZE);
    int v = 1;
    for (int z = 0; z < VSIZE; z++) {
        for (int y = 0; y < VSIZE; y++) {
            for (int x = 0; x < VSIZE; x++) {
                iprior_data1(1, v) = VAL * 1.5;
                iprior_data2(1, v) = VAL * 2.5;
                v++;
            }
        }
    }

    FabberRunData rundata;
    rundata.SetVoxelData("PSP_byname1_image", iprior_data1);
    rundata.SetVoxelData("PSP_byname2_image", iprior_data2);
    rundata.Set("param-spatial-priors", "II");

    PriorType prior(0, "a", rundata);
    ASSERT_EQ("a", prior.m_param_name);
    ASSERT_EQ(0, prior.m_idx);
    ASSERT_EQ('I', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("PSP_byname1_image", prior.m_filename);
    ASSERT_EQ(VSIZE * VSIZE * VSIZE, prior.m_image.Ncols());

    prior = PriorType(1, "z", rundata);
    ASSERT_EQ("z", prior.m_param_name);
    ASSERT_EQ(1, prior.m_idx);
    ASSERT_EQ('I', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("PSP_byname2_image", prior.m_filename);
    ASSERT_EQ(VSIZE * VSIZE * VSIZE, prior.m_image.Ncols());
}

TEST_F(PriorsTest, ImagePriorByName)
{
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix iprior_data1, iprior_data2;
    iprior_data1.ReSize(1, VSIZE * VSIZE * VSIZE);
    iprior_data2.ReSize(1, VSIZE * VSIZE * VSIZE);
    int v = 1;
    for (int z = 0; z < VSIZE; z++) {
        for (int y = 0; y < VSIZE; y++) {
            for (int x = 0; x < VSIZE; x++) {
                iprior_data1(1, v) = VAL * 1.5;
                iprior_data2(1, v) = VAL * 2.5;
                v++;
            }
        }
    }

    FabberRunData rundata;
    rundata.SetVoxelData("PSP_byname1_image", iprior_data1);
    rundata.Set("PSP_byname1", "a");
    rundata.Set("PSP_byname1_type", "I");

    PriorType prior(0, "a", rundata);
    ASSERT_EQ("a", prior.m_param_name);
    ASSERT_EQ(0, prior.m_idx);
    ASSERT_EQ('I', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("PSP_byname1_image", prior.m_filename);
    NEWMAT::RowVector image = prior.m_image;
    ASSERT_EQ(VSIZE * VSIZE * VSIZE, image.Ncols());
}

TEST_F(PriorsTest, ImagePriorByNameCorrectParam)
{
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix iprior_data1, iprior_data2;
    iprior_data1.ReSize(1, VSIZE * VSIZE * VSIZE);
    iprior_data2.ReSize(1, VSIZE * VSIZE * VSIZE);
    int v = 1;
    for (int z = 0; z < VSIZE; z++) {
        for (int y = 0; y < VSIZE; y++) {
            for (int x = 0; x < VSIZE; x++) {
                iprior_data1(1, v) = VAL * 1.5;
                iprior_data2(1, v) = VAL * 2.5;
                v++;
            }
        }
    }

    FabberRunData rundata;
    rundata.SetVoxelData("PSP_byname1_image", iprior_data1);
    rundata.Set("PSP_byname1", "b");
    rundata.Set("PSP_byname1_type", "I");

    PriorType prior(0, "a", rundata);
    ASSERT_EQ("a", prior.m_param_name);
    ASSERT_EQ(0, prior.m_idx);
    ASSERT_EQ('-', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("", prior.m_filename);
    ASSERT_EQ(0, prior.m_image.Ncols());

    prior = PriorType(1, "b", rundata);
    ASSERT_EQ("b", prior.m_param_name);
    ASSERT_EQ(1, prior.m_idx);
    ASSERT_EQ('I', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("PSP_byname1_image", prior.m_filename);
    ASSERT_EQ(VSIZE * VSIZE * VSIZE, prior.m_image.Ncols());
}

TEST_F(PriorsTest, ImagePriorByNameMultiple)
{
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix iprior_data1, iprior_data2;
    iprior_data1.ReSize(1, VSIZE * VSIZE * VSIZE);
    iprior_data2.ReSize(1, VSIZE * VSIZE * VSIZE);
    int v = 1;
    for (int z = 0; z < VSIZE; z++) {
        for (int y = 0; y < VSIZE; y++) {
            for (int x = 0; x < VSIZE; x++) {
                iprior_data1(1, v) = VAL * 1.5;
                iprior_data2(1, v) = VAL * 2.5;
                v++;
            }
        }
    }

    FabberRunData rundata;
    rundata.SetVoxelData("PSP_byname1_image", iprior_data1);
    rundata.SetVoxelData("PSP_byname2_image", iprior_data2);
    rundata.Set("PSP_byname1", "z");
    rundata.Set("PSP_byname1_type", "I");
    rundata.Set("PSP_byname2", "a");
    rundata.Set("PSP_byname2_type", "I");

    PriorType prior(0, "a", rundata);
    ASSERT_EQ("a", prior.m_param_name);
    ASSERT_EQ(0, prior.m_idx);
    ASSERT_EQ('I', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("PSP_byname2_image", prior.m_filename);
    ASSERT_EQ(VSIZE * VSIZE * VSIZE, prior.m_image.Ncols());

    prior = PriorType(1, "b", rundata);
    ASSERT_EQ("b", prior.m_param_name);
    ASSERT_EQ(1, prior.m_idx);
    ASSERT_EQ('-', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("", prior.m_filename);
    ASSERT_EQ(0, prior.m_image.Ncols());

    prior = PriorType(2, "z", rundata);
    ASSERT_EQ("z", prior.m_param_name);
    ASSERT_EQ(2, prior.m_idx);
    ASSERT_EQ('I', prior.m_type);
    ASSERT_EQ(-1, prior.m_prec);
    ASSERT_EQ("PSP_byname1_image", prior.m_filename);
    ASSERT_EQ(VSIZE * VSIZE * VSIZE, prior.m_image.Ncols());
}

TEST_F(PriorsTest, ImagePriorPrec)
{
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix iprior_data1, iprior_data2;
    iprior_data1.ReSize(1, VSIZE * VSIZE * VSIZE);
    iprior_data2.ReSize(1, VSIZE * VSIZE * VSIZE);
    int v = 1;
    for (int z = 0; z < VSIZE; z++) {
        for (int y = 0; y < VSIZE; y++) {
            for (int x = 0; x < VSIZE; x++) {
                iprior_data1(1, v) = VAL * 1.5;
                iprior_data2(1, v) = VAL * 2.5;
                v++;
            }
        }
    }

    FabberRunData rundata;
    rundata.SetVoxelData("PSP_byname1_image", iprior_data1);
    rundata.SetVoxelData("PSP_byname2_image", iprior_data2);
    rundata.Set("PSP_byname1", "z");
    rundata.Set("PSP_byname1_type", "I");
    rundata.Set("PSP_byname1_prec", "42");
    rundata.Set("PSP_byname2", "a");
    rundata.Set("PSP_byname2_type", "I");
    rundata.Set("PSP_byname2_prec", "24");

    PriorType prior(0, "a", rundata);
    ASSERT_EQ("a", prior.m_param_name);
    ASSERT_EQ(0, prior.m_idx);
    ASSERT_EQ('I', prior.m_type);
    ASSERT_FLOAT_EQ(24, prior.m_prec);
    ASSERT_EQ("PSP_byname2_image", prior.m_filename);
    ASSERT_EQ(VSIZE * VSIZE * VSIZE, prior.m_image.Ncols());

    prior = PriorType(2, "z", rundata);
    ASSERT_EQ("z", prior.m_param_name);
    ASSERT_EQ(2, prior.m_idx);
    ASSERT_EQ('I', prior.m_type);
    ASSERT_FLOAT_EQ(42, prior.m_prec);
    ASSERT_EQ("PSP_byname1_image", prior.m_filename);
    ASSERT_EQ(VSIZE * VSIZE * VSIZE, prior.m_image.Ncols());
}
