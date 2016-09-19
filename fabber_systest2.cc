#include <memory>
#include "newimage/newimageall.h"

#include "dataset.h"
#include "setup.h"
#include "fwdmodel.h"
#include "inference.h"

using namespace std;
using namespace NEWIMAGE;

#define X 5
#define Y 5
#define Z 5
#define NTIMES 106

// This test is a minimal test of the spatial VB method
int main()
{
    EasyLog::StartLog(".", true);

    // Create coordinates and data matrices
    Matrix voxelCoords, data;
    data.ReSize(NTIMES, X * Y * Z);
    voxelCoords.ReSize(3, X * Y * Z);
    int v = 1;
    for (int z = 0; z < Z; z++) {
        for (int y = 0; y < Y; y++) {
            for (int x = 0; x < X; x++) {
                voxelCoords(1, v) = x + 1;
                voxelCoords(2, v) = y + 1;
                voxelCoords(3, v) = z + 1;
                for (int n = 0; n < NTIMES; n++) {
                    data(n + 1, v) = 20 - ((n - (NTIMES / 2)) * (n - (NTIMES / 2)) / 400 + (x - (X / 2))
                            * (x - (X / 2)) + (y - (Y / 2)) * (y - (Y / 2)) + (z - (Z / 2)) * (z - (Z / 2)));
                }
                v++;
            }
        }
    }

    volume4D<float> data_out(X, Y, Z, 1);
    data_out.setmatrix(data);
    data_out.setDisplayMaximumMinimum(data_out.max(), data_out.min());
    save_volume4D(data_out, "data");

    // Run Fabber
    FabberRunData rundata;
    rundata.SetVoxelCoords(voxelCoords);
    rundata.SetMainVoxelData(data);
    rundata.Set("noise", "white");
    rundata.Set("model", "trivial");
    rundata.Set("method", "spatialvb");

    // Now run it:
    rundata.Run();

    // Lets get some data back and write it to a file
    Matrix mean = rundata.GetVoxelData("mean_p");
    volume4D<float> output(X, Y, Z, 1);
    output.setmatrix(mean);
    output.setDisplayMaximumMinimum(output.max(), output.min());
    save_volume4D(output, "mean_p");

    LOG << "fabber_library completed." << endl;

    return 0;
}
