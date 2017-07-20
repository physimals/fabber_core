//
// Tests specifically for the VB method
//
// Each test runs in both spatialvb and vb modes

#ifndef NO_NEWIMAGE
#include "newimage/newimageall.h"
#endif

#include "gtest/gtest.h"

#include "easylog.h"
#include "inference.h"
#include "inference_myvb.h"
#include "rundata_newimage.h"
#include "setup.h"

namespace
{

class VbTest : public ::testing::TestWithParam<string>
{
protected:
    VbTest()
        : rundata(NULL)
        , model(NULL)
        , vb(NULL)
    {
        FabberSetup::SetupDefaults();
    }

    virtual ~VbTest()
    {
        FabberSetup::Destroy();
    }

    virtual void SetUp()
    {
        rundata = new FabberRunDataNewimage();
        rundata->SetLogger(&log);
        rundata->Set("method", GetParam());
    }

    virtual void TearDown()
    {
        delete vb;
        delete rundata;
    }

    void Run()
    {
        std::auto_ptr<FwdModel> fwd_model(FwdModel::NewFromName(rundata->GetString("model")));
        fwd_model->Initialize(*rundata);
        std::auto_ptr<InferenceTechnique> infer(InferenceTechnique::NewFromName(rundata->GetString("method")));
        
        infer->Initialize(fwd_model.get(), *rundata);  
        infer->DoCalculations(*rundata);
        infer->SaveResults(*rundata);
    }

    EasyLog log;
    NEWMAT::Matrix voxelCoords;
    FabberRunData *rundata;
    FwdModel *model;
    Vb *vb;
};

// Test image priors. Note that this just
// checks the code works when they are specified
// not that they are actually having any effect!
TEST_P(VbTest, ImagePriors)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data, iprior_data;
    data.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    iprior_data.ReSize(1, VSIZE * VSIZE * VSIZE);
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
                iprior_data(1, v) = VAL * 1.5;
                v++;
            }
        }
    }
    rundata->SetVoxelCoords(voxelCoords);
    rundata->SetVoxelData("data", data);
    rundata->Set("noise", "white");
    rundata->Set("model", "poly");
    rundata->Set("degree", "0");

    rundata->Set("PSP_byname1", "c0");
    rundata->Set("PSP_byname1_type", "I");
    rundata->SetVoxelData("PSP_byname1_image", iprior_data);
    Run();
}

// Test image priors with a precision which is to
// high for convergence
TEST_P(VbTest, ImagePriorsPrecTooHigh)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data, iprior_data;
    data.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    iprior_data.ReSize(1, VSIZE * VSIZE * VSIZE);
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
                iprior_data(1, v) = VAL * 1.5;
                v++;
            }
        }
    }
    rundata->SetVoxelCoords(voxelCoords);
    rundata->SetVoxelData("data", data);
    rundata->Set("noise", "white");
    rundata->Set("model", "poly");
    rundata->Set("degree", "0");
    rundata->Set("method", GetParam());

    rundata->Set("PSP_byname1", "c0");
    rundata->Set("PSP_byname1_type", "I");
    rundata->SetVoxelData("PSP_byname1_image", iprior_data);
    rundata->Set("PSP_byname1_prec", "1e12");
    Run();

    NEWMAT::Matrix mean = rundata->GetVoxelData("mean_c0");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        // we specified an image prior with high precision
        // so we expect the output to be close to that
        ASSERT_NEAR(mean(1, i + 1), VAL * 1.5, VAL*0.1);
    }
}

// Test image priors with a precision that is low enough
// to get to the 'right' answer
TEST_P(VbTest, ImagePriorsPrecLow)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data, iprior_data;
    data.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    iprior_data.ReSize(1, VSIZE * VSIZE * VSIZE);
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
                iprior_data(1, v) = VAL * 1.5;
                v++;
            }
        }
    }
    rundata->SetVoxelCoords(voxelCoords);
    rundata->SetVoxelData("data", data);
    rundata->Set("noise", "white");
    rundata->Set("model", "poly");
    rundata->Set("degree", "0");
    rundata->Set("method", GetParam());

    rundata->Set("PSP_byname1", "c0");
    rundata->Set("PSP_byname1_type", "I");
    rundata->SetVoxelData("PSP_byname1_image", iprior_data);
    rundata->Set("PSP_byname1_prec", "1e-5");
    Run();

    NEWMAT::Matrix mean = rundata->GetVoxelData("mean_c0");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        // We expect the parameter to be about 'right' because
        // we specified an image prior with low precision
        ASSERT_NEAR(mean(1, i + 1), VAL * 2, VAL*0.1);
    }
}

#ifndef NO_NEWIMAGE
// Test image priors when stored in a file
TEST_P(VbTest, ImagePriorsFile)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;
    string FILENAME = "imageprior_data_temp";

    // Create coordinates and data matrices
    NEWMAT::Matrix voxelCoords, data, iprior_data;
    data.ReSize(NTIMES, VSIZE * VSIZE * VSIZE);
    iprior_data.ReSize(1, VSIZE * VSIZE * VSIZE);
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
                iprior_data(1, v) = VAL * 1.5;
                v++;
            }
        }
    }

    // Save image prior data to a file
    NEWIMAGE::volume4D<float> data_out(VSIZE, VSIZE, VSIZE, 1);
    data_out.setmatrix(iprior_data);
    data_out.setDisplayMaximumMinimum(data_out.max(), data_out.min());
    save_volume4D(data_out, FILENAME);
    iprior_data.ReSize(1, 1);

    rundata->SetVoxelCoords(voxelCoords);
    rundata->SetVoxelData("data", data);
    rundata->Set("noise", "white");
    rundata->Set("model", "poly");
    rundata->Set("degree", "0");
    rundata->Set("method", GetParam());

    rundata->Set("PSP_byname1", "c0");
    rundata->Set("PSP_byname1_type", "I");
    rundata->Set("PSP_byname1_image", FILENAME);
    rundata->Set("PSP_byname1_prec", "1e12");

    Run();

    NEWMAT::Matrix mean = rundata->GetVoxelData("mean_c0");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        ASSERT_NEAR(mean(1, i + 1), VAL * 1.5, VAL*0.1);
    }

    remove(string(FILENAME + ".nii.gz").c_str());
}
#endif

// Test restarting VB run
TEST_P(VbTest, Restart)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;
    int REPEATS = 50;
    int DEGREE = 5;

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

    // Do just 1 iteration
    rundata->SetVoxelCoords(voxelCoords);
    rundata->SetVoxelData("data", data);
    rundata->Set("noise", "white");
    rundata->Set("model", "poly");
    rundata->Set("degree", stringify(DEGREE));
    rundata->Set("method", "vb");
    rundata->Set("max-iterations", "1");
    Run();

    // Make sure not converged after first iteration!
    NEWMAT::Matrix mean = rundata->GetVoxelData("mean_c0");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    ASSERT_NE(VAL, mean(1, 1));

    mean = rundata->GetVoxelData("mean_c2");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    ASSERT_NE(VAL * 1.5, mean(1, 1));

    // Save results and repeat. Should converge even though we
    // only do 1 iteration each time.
    for (int repeat = 0; repeat < REPEATS; repeat++)
    {
        NEWMAT::Matrix mvns = rundata->GetVoxelData("finalMVN");
        //ASSERT_EQ(mvns.Nrows(), 7);
        // This was just so you could see the convergence
        //mean = rundata->GetVoxelData("mean_c0");
        //cout << mean(1, 1) << " != " << VAL << endl;
        //mean = rundata->GetVoxelData("mean_c2");
        //cout << mean(1, 1) << " != " << VAL*1.5 << endl;

        TearDown();
        SetUp();

        rundata->Set("noise", "white");
        rundata->Set("model", "poly");
        rundata->Set("degree", stringify(DEGREE));
        rundata->Set("method", "vb");
        rundata->Set("max-iterations", "1");
        rundata->Set("continue-from-mvn", "mvns");
        rundata->SetVoxelCoords(voxelCoords);
        rundata->SetVoxelData("data", data);
        rundata->SetVoxelData("mvns", mvns);

        Run();
    }

    mean = rundata->GetVoxelData("mean_c0");

    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        ASSERT_FLOAT_EQ(VAL, mean(1, i + 1));
    }
    mean = rundata->GetVoxelData("mean_c1");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        // GTEST has difficulty with comparing floats to 0
        ASSERT_FLOAT_EQ(1, mean(1, i + 1) + 1);
    }
    mean = rundata->GetVoxelData("mean_c2");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        ASSERT_FLOAT_EQ(VAL * 1.5, mean(1, i + 1));
    }
}

// Test restarting VB run with the output-only option
TEST_P(VbTest, RestartOutputOnly)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;
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

    rundata->SetVoxelCoords(voxelCoords);
    rundata->SetVoxelData("data", data);
    rundata->Set("noise", "white");
    rundata->Set("model", "poly");
    rundata->Set("degree", stringify(DEGREE));
    rundata->Set("method", "vb");
    rundata->Set("max-iterations", "100");
    rundata->SetBool("save-mean", false);
    rundata->SetBool("save-modelfit", false);
    rundata->SetBool("save-mvn", true);
    Run();

    NEWMAT::Matrix mvns = rundata->GetVoxelData("finalMVN");
    TearDown();
    SetUp();
    rundata->Set("noise", "white");
    rundata->Set("model", "poly");
    rundata->Set("degree", stringify(DEGREE));
    rundata->Set("method", "vb");
    rundata->SetBool("save-mean", true);
    rundata->SetBool("save-model-fit", true);
    rundata->SetBool("save-mvn", false);
    rundata->SetBool("output-only", true);
    rundata->Set("continue-from-mvn", "mvns");
    rundata->SetVoxelCoords(voxelCoords);
    rundata->SetVoxelData("data", data);
    rundata->SetVoxelData("mvns", mvns);
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

#ifndef NO_NEWIMAGE
// Test restarting VB run from a file
TEST_P(VbTest, RestartFromFile)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;
    int REPEATS = 50;
    int DEGREE = 5;
    string FILENAME = "temp_mvns";

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

    // Do just 1 iteration
    rundata->SetVoxelCoords(voxelCoords);
    rundata->SetVoxelData("data", data);
    rundata->Set("noise", "white");
    rundata->Set("model", "poly");
    rundata->Set("degree", stringify(DEGREE));
    rundata->Set("method", "vb");
    rundata->Set("max-iterations", "1");
    Run();

    // Make sure not converged after first iteration!
    NEWMAT::Matrix mean = rundata->GetVoxelData("mean_c0");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    ASSERT_NE(VAL, mean(1, 1));

    mean = rundata->GetVoxelData("mean_c2");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    ASSERT_NE(VAL * 1.5, mean(1, 1));

    // Save results and repeat. Should converge even though we
    // only do 1 iteration each time.
    for (int repeat = 0; repeat < REPEATS; repeat++)
    {
        NEWMAT::Matrix mvns = rundata->GetVoxelData("finalMVN");
        NEWIMAGE::volume4D<float> data_out(VSIZE, VSIZE, VSIZE, mvns.Nrows());
        data_out.setmatrix(mvns);
        data_out.setDisplayMaximumMinimum(data_out.max(), data_out.min());
        save_volume4D(data_out, FILENAME);
        mvns.ReSize(1, 1);

        TearDown();
        SetUp();

        //ASSERT_EQ(mvns.Nrows(), 7);
        rundata->SetVoxelCoords(voxelCoords);
        rundata->SetVoxelData("data", data);
        rundata->Set("max-iterations", "1");
        rundata->Set("continue-from-mvn", FILENAME);
        rundata->Set("noise", "white");
        rundata->Set("model", "poly");
        rundata->Set("degree", stringify(DEGREE));
        rundata->Set("method", "vb");
        // This was just so you could see the convergence
        //mean = rundata->GetVoxelData("mean_c0");
        //cout << mean(1, 1) << " != " << VAL << endl;
        //mean = rundata->GetVoxelData("mean_c2");
        //cout << mean(1, 1) << " != " << VAL*1.5 << endl;

        Run();

        // Don't pick up last iteration run
        remove((FILENAME + ".nii.gz").c_str());
    }

    mean = rundata->GetVoxelData("mean_c0");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        ASSERT_FLOAT_EQ(VAL, mean(1, i + 1));
    }
    mean = rundata->GetVoxelData("mean_c1");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        // GTEST has difficulty with comparing floats to 0
        ASSERT_FLOAT_EQ(1, mean(1, i + 1) + 1);
    }
    mean = rundata->GetVoxelData("mean_c2");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        ASSERT_FLOAT_EQ(VAL * 1.5, mean(1, i + 1));
    }
}
#endif

// Test fitting to a simple polynomial model with noise
TEST_P(VbTest, ArNoise)
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
                    float noise = (float(rand()) / RAND_MAX - 0.5) * VAL / 200;
                    data(n + 1, v) = VAL + (1.5 * VAL) * (n + 1) * (n + 1) - 2 * VAL * (n + 1) * (n + 1) * (n + 1) + noise;
                }
                v++;
            }
        }
    }

    // Do just 1 iteration
    rundata->SetVoxelCoords(voxelCoords);
    rundata->SetVoxelData("data", data);
    rundata->Set("noise", "ar");
    rundata->Set("model", "poly");
    rundata->Set("max-iterations", "50");
    rundata->Set("degree", stringify(DEGREE));
    rundata->Set("method", GetParam());
    rundata->Run();

    NEWMAT::Matrix mean = rundata->GetVoxelData("mean_c0");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), n_voxels);
    for (int i = 0; i < n_voxels; i++)
    {
        EXPECT_NEAR(VAL, mean(1, i + 1), 0.2);
    }

    mean = rundata->GetVoxelData("mean_c1");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), n_voxels);
    for (int i = 0; i < n_voxels; i++)
    {
        // GTEST has difficulty with comparing floats to 0
        EXPECT_NEAR(0, mean(1, i + 1), 0.2);
    }

    mean = rundata->GetVoxelData("mean_c2");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), n_voxels);
    for (int i = 0; i < n_voxels; i++)
    {
        EXPECT_NEAR(VAL * 1.5, mean(1, i + 1), 0.2);
    }

    mean = rundata->GetVoxelData("mean_c3");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), n_voxels);
    for (int i = 0; i < n_voxels; i++)
    {
        EXPECT_NEAR(-VAL * 2, mean(1, i + 1), 0.2);
    }
}

// Test fitting to a simple polynomial model with noise
TEST_P(VbTest, WhiteNoise)
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
                    float noise = (float(rand()) / RAND_MAX - 0.5) * VAL / 100;
                    data(n + 1, v) = VAL + (1.5 * VAL) * (n + 1) * (n + 1) - 2 * VAL * (n + 1) * (n + 1) * (n + 1) + noise;
                }
                v++;
            }
        }
    }

    // Do just 1 iteration
    rundata->SetVoxelCoords(voxelCoords);
    rundata->SetVoxelData("data", data);
    rundata->Set("noise", "white");
    rundata->Set("model", "poly");
    rundata->Set("max-iterations", "50");
    rundata->Set("degree", stringify(DEGREE));
    rundata->Set("method", GetParam());
    rundata->Run();

    NEWMAT::Matrix mean = rundata->GetVoxelData("mean_c0");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), n_voxels);
    for (int i = 0; i < n_voxels; i++)
    {
        EXPECT_NEAR(VAL, mean(1, i + 1), 0.2);
    }

    mean = rundata->GetVoxelData("mean_c1");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), n_voxels);
    for (int i = 0; i < n_voxels; i++)
    {
        // GTEST has difficulty with comparing floats to 0
        EXPECT_NEAR(0, mean(1, i + 1), 0.2);
    }

    mean = rundata->GetVoxelData("mean_c2");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), n_voxels);
    for (int i = 0; i < n_voxels; i++)
    {
        EXPECT_NEAR(VAL * 1.5, mean(1, i + 1), 0.2);
    }

    mean = rundata->GetVoxelData("mean_c3");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), n_voxels);
    for (int i = 0; i < n_voxels; i++)
    {
        EXPECT_NEAR(-VAL * 2, mean(1, i + 1), 0.2);
    }
}

#ifdef __FABBER_MOTION
// Turn motion correction on, but no motion to correct!
TEST_P(VbTest, MotionCorNull)
{
    int NTIMES = 10; // needs to be even
    int VSIZE = 5;
    float VAL = 7.32;
    int REPEATS = 50;
    int DEGREE = 5;
    string FILENAME = "temp_mvns";

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

    // Do just 1 iteration
    rundata->SetVoxelCoords(voxelCoords);
    rundata->SetVoxelData("data", data);
    rundata->Set("noise", "white");
    rundata->Set("mcsteps", "5");
    rundata->Set("model", "poly");
    rundata->Set("degree", stringify(DEGREE));
    rundata->Set("method", "vb");
    rundata->Set("max-iterations", "50");
    Run();

    NEWMAT::Matrix mean = rundata->GetVoxelData("mean_c0");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        ASSERT_FLOAT_EQ(VAL, mean(1, i + 1));
    }
    mean = rundata->GetVoxelData("mean_c1");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        // GTEST has difficulty with comparing floats to 0
        ASSERT_FLOAT_EQ(1, mean(1, i + 1) + 1);
    }
    mean = rundata->GetVoxelData("mean_c2");
    ASSERT_EQ(mean.Nrows(), 1);
    ASSERT_EQ(mean.Ncols(), VSIZE * VSIZE * VSIZE);
    for (int i = 0; i < VSIZE * VSIZE * VSIZE; i++)
    {
        ASSERT_FLOAT_EQ(VAL * 1.5, mean(1, i + 1));
    }
}
#endif

INSTANTIATE_TEST_CASE_P(VbTests, VbTest, ::testing::Values("vb", "spatialvb"));
}
