#include "gtest/gtest.h"

#include "priors.h"

#include "rundata.h"

// Tests for the base class - static ExpandPriorTypesString is the only relevant method

class PriorTest : public ::testing::Test
{
};

TEST_F(PriorTest, NoPriorsSpecified)
{
    unsigned int NUM_PARAMS = 4;
    string prior_types = Prior::ExpandPriorTypesString("", NUM_PARAMS);
    ASSERT_EQ(NUM_PARAMS, prior_types.size());
    for (unsigned int i = 0; i < NUM_PARAMS; i++)
    {
        ASSERT_EQ('-', prior_types[i]);
    }
}

TEST_F(PriorTest, AllPriorsSpecified)
{
    string TYPES_STRING = "ABCD";

    string prior_types = Prior::ExpandPriorTypesString(TYPES_STRING, TYPES_STRING.size());
    ASSERT_EQ(TYPES_STRING, prior_types);
}

TEST_F(PriorTest, SomePriorsSpecified)
{
    string TYPES_STRING = "ABCD";
    unsigned int NUM_EXTRA = 3;

    string prior_types
        = Prior::ExpandPriorTypesString(TYPES_STRING, TYPES_STRING.size() + NUM_EXTRA);
    for (unsigned int i = 0; i < TYPES_STRING.size(); i++)
    {
        ASSERT_EQ(TYPES_STRING[i], prior_types[i]);
    }
    // Extras should be model default
    for (unsigned int i = 0; i < NUM_EXTRA; i++)
    {
        ASSERT_EQ('-', prior_types[TYPES_STRING.size() + i]);
    }
}

TEST_F(PriorTest, TrailingPlus)
{
    string TYPES_STRING = "ABC+";
    unsigned int NUM_PARAMS = 7;

    string prior_types = Prior::ExpandPriorTypesString(TYPES_STRING, NUM_PARAMS);
    ASSERT_EQ(prior_types.size(), NUM_PARAMS);
    for (unsigned int i = 0; i < NUM_PARAMS; i++)
    {
        if (i < TYPES_STRING.size() - 1)
            ASSERT_EQ(TYPES_STRING[i], prior_types[i]);
        else
            ASSERT_EQ(TYPES_STRING[TYPES_STRING.size() - 2], prior_types[i]);
    }
}

TEST_F(PriorTest, EmbeddedPlus)
{
    string TYPES_STRING = "AB+C";
    unsigned int NUM_PARAMS = 7;

    string prior_types = Prior::ExpandPriorTypesString(TYPES_STRING, NUM_PARAMS);
    ASSERT_EQ(prior_types.size(), NUM_PARAMS);
    ASSERT_EQ(prior_types, "ABBBBBC");
}

TEST_F(PriorTest, PointlessTrailingPlus)
{
    string TYPES_STRING = "ABCDE+";
    unsigned int NUM_PARAMS = TYPES_STRING.size() - 1;

    string prior_types = Prior::ExpandPriorTypesString(TYPES_STRING, NUM_PARAMS);
    ASSERT_EQ(prior_types.size(), NUM_PARAMS);
    for (unsigned int i = 0; i < NUM_PARAMS; i++)
    {
        ASSERT_EQ(TYPES_STRING[i], prior_types[i]);
    }
}

TEST_F(PriorTest, LeadingPlus)
{
    string TYPES_STRING = "+ABCDE";
    unsigned int NUM_PARAMS = TYPES_STRING.size() + 3;

    string prior_types = Prior::ExpandPriorTypesString(TYPES_STRING, NUM_PARAMS);
    ASSERT_EQ(prior_types.size(), NUM_PARAMS);
    ASSERT_EQ(prior_types, "----ABCDE");
}

TEST_F(PriorTest, PointlessLeadingPlus)
{
    string TYPES_STRING = "+ABCDE";
    unsigned int NUM_PARAMS = TYPES_STRING.size() - 1;

    string prior_types = Prior::ExpandPriorTypesString(TYPES_STRING, NUM_PARAMS);
    ASSERT_EQ(prior_types.size(), NUM_PARAMS);
    for (unsigned int i = 0; i < NUM_PARAMS; i++)
    {
        ASSERT_EQ(TYPES_STRING[i + 1], prior_types[i]);
    }
}

TEST_F(PriorTest, PlusOnly)
{
    string TYPES_STRING = "+";
    unsigned int NUM_PARAMS = 12;

    string prior_types = Prior::ExpandPriorTypesString(TYPES_STRING, NUM_PARAMS);
    ASSERT_EQ(prior_types.size(), NUM_PARAMS);
    for (unsigned int i = 0; i < NUM_PARAMS; i++)
    {
        ASSERT_EQ('-', prior_types[i]);
    }
}

TEST_F(PriorTest, TooManyPriorsSpecified)
{
    string TYPES_STRING = "ABCD";
    ASSERT_THROW(
        Prior::ExpandPriorTypesString(TYPES_STRING, TYPES_STRING.size() - 1), std::exception);
}

TEST_F(PriorTest, TooManyPlusses)
{
    string TYPES_STRING = "A+BCD+";
    ASSERT_THROW(Prior::ExpandPriorTypesString(TYPES_STRING, TYPES_STRING.size()), std::exception);
}

TEST_F(PriorTest, TooManyPlussesEnd)
{
    string TYPES_STRING = "ABCD++";
    ASSERT_THROW(Prior::ExpandPriorTypesString(TYPES_STRING, TYPES_STRING.size()), std::exception);
}

// Tests for concrete subclasses

#define PARAM_NAME "SomeParameter"
#define PARAM_IDX 3
#define PRIOR_MEAN 3.9
#define PRIOR_VAR 5.3
#define POST_MEAN 6.3
#define POST_VAR 2.3
#define PRIOR_TYPE 'Q'
#define NUM_VOXELS 17
#define IMAGE_FNAME "IMAGE"

class DefaultPriorTest : public ::testing::Test
{
};

TEST_F(DefaultPriorTest, BasicProperties)
{
    Parameter p(PARAM_IDX, PARAM_NAME, DistParams(PRIOR_MEAN, PRIOR_VAR),
        DistParams(POST_MEAN, POST_VAR), PRIOR_TYPE);
    DefaultPrior prior(p);
    ASSERT_EQ(prior.m_param_name, PARAM_NAME);
    ASSERT_EQ(prior.m_idx, PARAM_IDX);
    ASSERT_EQ(prior.m_type_code, PRIOR_TYPE);
    ASSERT_EQ(prior.m_params.mean(), PRIOR_MEAN);
    ASSERT_EQ(prior.m_params.var(), PRIOR_VAR);
}

TEST_F(DefaultPriorTest, IgnoresTransform)
{
    Parameter p(PARAM_IDX, PARAM_NAME, DistParams(PRIOR_MEAN, PRIOR_VAR),
        DistParams(POST_MEAN, POST_VAR), PRIOR_TYPE, TRANSFORM_LOG());
    DefaultPrior prior(p);
    ASSERT_EQ(prior.m_param_name, PARAM_NAME);
    ASSERT_EQ(prior.m_idx, PARAM_IDX);
    ASSERT_EQ(prior.m_type_code, PRIOR_TYPE);
    ASSERT_EQ(prior.m_params.mean(), PRIOR_MEAN);
    ASSERT_EQ(prior.m_params.var(), PRIOR_VAR);
}

TEST_F(DefaultPriorTest, ApplyToMVN)
{
    Parameter p(PARAM_IDX, PARAM_NAME, DistParams(PRIOR_MEAN, PRIOR_VAR),
        DistParams(POST_MEAN, POST_VAR), PRIOR_TYPE);
    DefaultPrior prior(p);

    MVNDist mvn(PARAM_IDX + 7);
    RunContext ctx(1);
    prior.ApplyToMVN(&mvn, ctx);
    for (int i = 1; i <= NUM_VOXELS; i++)
    {
        ctx.v = i;
        prior.ApplyToMVN(&mvn, ctx);
        ASSERT_EQ(mvn.means(PARAM_IDX + 1), PRIOR_MEAN);
        ASSERT_EQ(mvn.GetCovariance()(PARAM_IDX + 1, PARAM_IDX + 1), PRIOR_VAR);
    }
}

class ImagePriorTest : public ::testing::Test
{
};

TEST_F(ImagePriorTest, BasicProperties)
{
    NEWMAT::Matrix data;
    data.ReSize(1, NUM_VOXELS);
    for (int i = 1; i <= NUM_VOXELS; i++)
        data(1, i) = i * 2.3 - 7.4;

    FabberRunData rundata;
    rundata.SetVoxelData(IMAGE_FNAME, data);

    Parameter p(PARAM_IDX, PARAM_NAME, DistParams(PRIOR_MEAN, PRIOR_VAR),
        DistParams(POST_MEAN, POST_VAR), PRIOR_TYPE);
    p.options["image"] = IMAGE_FNAME;

    ImagePrior prior(p, rundata);

    ASSERT_EQ(prior.m_param_name, PARAM_NAME);
    ASSERT_EQ(prior.m_idx, PARAM_IDX);
    ASSERT_EQ(prior.m_type_code, PRIOR_TYPE);
    ASSERT_EQ(prior.m_params.mean(), PRIOR_MEAN);
    ASSERT_EQ(prior.m_params.var(), PRIOR_VAR);
}

TEST_F(ImagePriorTest, ApplyToMVN)
{
    NEWMAT::Matrix data;
    data.ReSize(1, NUM_VOXELS);
    for (int i = 1; i <= NUM_VOXELS; i++)
        data(1, i) = i * 2.3 - 7.4;

    FabberRunData rundata;
    rundata.SetVoxelData(IMAGE_FNAME, data);

    Parameter p(PARAM_IDX, PARAM_NAME, DistParams(PRIOR_MEAN, PRIOR_VAR),
        DistParams(POST_MEAN, POST_VAR), PRIOR_TYPE);
    p.options["image"] = IMAGE_FNAME;
    ImagePrior prior(p, rundata);

    MVNDist mvn(PARAM_IDX + 7);
    RunContext ctx(NUM_VOXELS);
    for (int i = 1; i <= NUM_VOXELS; i++)
    {
        ctx.v = i;
        prior.ApplyToMVN(&mvn, ctx);
        ASSERT_EQ(mvn.means(PARAM_IDX + 1), data(1, i));
        ASSERT_EQ(mvn.GetCovariance()(PARAM_IDX + 1, PARAM_IDX + 1), PRIOR_VAR);
    }
}
