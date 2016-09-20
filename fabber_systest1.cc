#include <memory>
#include "newimage/newimageall.h"

#include "dataset.h"
#include "setup.h"
#include "fwdmodel.h"
#include "inference.h"

using namespace std;
using namespace NEWIMAGE;

/**
 * Convert a mask specifying voxels of interest into a list of xyz co-ordinates of
 * voxels of interest
 *
 * @param mask Input mask volume whose entries are 1 for voxels of interest, 0 otherwise
 * @param voxelCoords Return matrix in which every column is a the xyz vector or a 
 *                    voxel's co-ordinates
 */
void ConvertMaskToVoxelCoordinates(const volume<float>& mask, Matrix& voxelCoords)
{
	Tracer_Plus("ConvertMaskToVoxelCoordinates");

	// Mask has previously been binarized to 0 or 1 to identify
	// voxels of interest, so sum is count of these voxels
	ColumnVector preThresh((int) mask.sum());
	const int nVoxels = preThresh.Nrows();

	// Populate preThresh with offset of voxel into matrix
	// starting at 0 FIXME repeated from calculateNeigbours
	int offset(0);
	int count(1);
	for (int z = 0; z < mask.zsize(); z++)
	{
		for (int y = 0; y < mask.ysize(); y++)
		{
			for (int x = 0; x < mask.xsize(); x++)
			{
				if (mask(x, y, z) != 0)
				{
					preThresh(count++) = offset;
				}
				offset++;
			}
		}
	}

	// Dimensions of mask matrix
	const int dims[3] =
	{ mask.xsize(), mask.ysize(), mask.zsize() };

	// Basic sanity checks. Note that xdim/ydim/zdim are the physical
	// dimensions of the voxels
	assert(mask.xsize()*mask.ysize()*mask.zsize() > 0);
	assert(mask.xdim()*mask.ydim()*mask.zdim() > 0);

	LOG << "Calculating distance matrix, using voxel dimensions: " << mask.xdim() << " by " << mask.ydim() << " by "
			<< mask.zdim() << " mm\n";

	ColumnVector positions[3]; // indices
	positions[0].ReSize(nVoxels);
	positions[1].ReSize(nVoxels);
	positions[2].ReSize(nVoxels);

	// Go through each voxel. Note that preThresh is a ColumnVector and indexes from 1 not 0
	for (int vox = 1; vox <= nVoxels; vox++)
	{
		int pos = (int) preThresh(vox) - 1; // preThresh appears to be 1-indexed. FIXME not sure, from above looks like 0 offset is possible for first voxel
		assert(pos>=0);

		positions[0](vox) = pos % dims[0]; // X co-ordinate of offset from preThresh
		pos = pos / dims[0];
		positions[1](vox) = pos % dims[1]; // Y co-ordinate of offset from preThresh
		pos = pos / dims[1];
		positions[2](vox) = pos % dims[2]; // Z co-ordinate of offset from preThresh
		pos = pos / dims[2];
		//LOG << vox << ' ' << (int)preThresh(vox)<< ' ' << dims[0] << ' ' << dims[1] << ' ' << dims[2] << ' ' << pos << endl;
		assert(pos == 0); // fails if preThresh >= product of dims[0,1,2]

		// FIXME looks like old code below - remove?

		//double pos = preThresh(vox);
		//positions[2](vox) = floor(pos/dims[0]/dims[1]);
		//pos -= positions[2](vox)*dims[0]*dims[1];
		//positions[1](vox) = floor(pos/dims[0]);
		//pos -= positions[1](vox)*dims[0];
		//positions[0](vox) = pos;
		//assert(positions[2](vox) < dims[2]);
		//assert(positions[1](vox) < dims[1]);
		//assert(positions[0](vox) < dims[0]);
		assert(preThresh(vox)-1 == positions[0](vox) + positions[1](vox)*dims[0] + positions[2](vox)*dims[0]*dims[1] );
	}

	// Turn 3 column vectors into a matrix where each column is a voxel. note that
	// need to transpose to do this
	voxelCoords = (positions[0] | positions[1] | positions[2]).t();
}

// A simple test program to show how to use the library
int main()
{
	EasyLog::StartLog(".", true);

	LOG << "FABBER_system test 1" << endl;

	FabberRunData rundata;

	LOG << "Loading mask from mask.nii.gz..." << endl;
	rundata.Set("mask", "mask");
	LOG << "Loading data from data.nii.gz..." << endl;
	rundata.Set("data", "data");

	rundata.Set("noise", (string) "white");
	rundata.Set("basis", (string) "design.mat");

	// Now run it:
	rundata.Run();

	// Lets get some data back and write it to a file
	// Could have done this automatically using
	// rundata.SetSaveData(true)

	Matrix mean = rundata.GetVoxelData("mean_Parameter_1");
	volume4D<float> output(mask.xsize(), mask.ysize(), mask.zsize(), 1);
	output.setmatrix(mean, mask);
	output.setDisplayMaximumMinimum(output.max(), output.min());
	save_volume4D(output, "mean_Parameter_1");

	LOG << "fabber_library completed." << endl;

	return 0;
}
