#include "fabber_io_newimage.h"

#include "dataset.h"
#include "easylog.h"

#include "newimage/newimage.h"
#include "newimage/newimageio.h"
#include "newmat.h"

#include <string>
#include <vector>
#include <ostream>

using namespace std;
using namespace NEWIMAGE;
using NEWMAT::Matrix;

static void DumpVolumeInfo(const volume4D<float>& info, ostream& out = LOG)
{
	LOG << "FabberIoNewimage::Dimensions: x=" << info.xsize() << ", y=" << info.ysize() << ", z=" << info.zsize()
			<< ", vols=" << info.tsize() << endl;
	LOG << "FabberIoNewimage::Voxel size: x=" << info.xdim() << "mm, y=" << info.ydim() << "mm, z=" << info.zdim()
			<< "mm, TR=" << info.tdim() << " sec\n";
	LOG << "FabberIoNewimage::Intents: " << info.intent_code() << ", " << info.intent_param(1) << ", "
			<< info.intent_param(2) << ", " << info.intent_param(3) << endl;
}

static void DumpVolumeInfo(const volume<float>& info, ostream& out = LOG)
{
	LOG << "FabberIoNewimage::Dimensions: x=" << info.xsize() << ", y=" << info.ysize() << ", z=" << info.zsize()
			<< ", vols=1" << endl;
	LOG << "FabberIoNewimage::Voxel size: x=" << info.xdim() << "mm, y=" << info.ydim() << "mm, z=" << info.zdim()
			<< "mm, TR=1" << " sec\n";
	LOG << "FabberIoNewimage::Intents: " << info.intent_code() << ", " << info.intent_param(1) << ", "
			<< info.intent_param(2) << ", " << info.intent_param(3) << endl;
}

FabberIoNewimage::FabberIoNewimage() :
		m_have_mask(false), m_have_coords(false)
{
}

void FabberIoNewimage::LoadMask(std::string filename)
{
	LOG_ERR("FabberIoNewimage::Loading mask data from '" + filename << "'" << endl);
	read_volume(m_mask, filename);
	m_mask.binarise(1e-16, m_mask.max() + 1, exclusive);
	DumpVolumeInfo(m_mask);
	m_have_mask = true;

	SetVoxelCoordsFromExtent(m_mask.xsize(), m_mask.ysize(), m_mask.zsize());
}

Matrix FabberIoNewimage::LoadVoxelData(std::string filename)
{
	// Load the data file using Newimage library
	LOG << "FabberIoNewimage::Loading data from '" + filename << "'" << endl;
	volume4D<float> vol;
	Matrix ret;
	try
	{
		read_volume4D(vol, filename);
	} catch (Exception &e)
	{
		throw DataLoadError(filename);
	}
	DumpVolumeInfo (vol);

	if (!m_have_coords)
		SetVoxelCoordsFromExtent(vol.xsize(), vol.ysize(), vol.zsize());

	try
	{
		if (m_have_mask)
		{
			LOG << "     Applying mask to data..." << endl;
			ret = vol.matrix(m_mask);
		}
		else
		{
			ret = vol.matrix();
		}
	} catch (exception &e)
	{
		LOG << "*** NEWMAT error while thresholding time-series... Most likely a dimension mismatch. ***\n";
		throw e;
	}
	return ret;
}

void FabberIoNewimage::SaveVoxelData(NEWMAT::Matrix &data, vector<int> extent, string filename, VoxelDataType data_type)
{
	LOG << "FabberIoNewimage::Saving to nifti: " << filename << endl;
	int nifti_intent_code;
	switch(data_type) {
	case VDT_MVN:
		nifti_intent_code=NIFTI_INTENT_SYMMATRIX;
		break;
	default:
		nifti_intent_code=NIFTI_INTENT_NONE;
	}

	int data_size = data.Nrows();
	volume4D<float> output(extent[0], extent[1], extent[2], data_size);
	if (m_have_mask)
	{
		output.setmatrix(data, m_mask);
	}
	else
	{
		output.setmatrix(data);
	}

	output.set_intent(nifti_intent_code, 0, 0, 0);
	output.setDisplayMaximumMinimum(output.max(), output.min());

	// FIXME need to use logger to get outdir as this does the ++ appending, however
	// this assumes use of CL tool and log to file
	string outDir = EasyLog::GetOutputDirectory();
	if (outDir != "")
		filename = outDir + "/" + filename;
	save_volume4D(output, filename);
}

Matrix FabberIoNewimage::GetVoxelCoords()
{
	return m_coords;
}

void FabberIoNewimage::SetVoxelCoordsFromExtent(int nx, int ny, int nz)
{
	LOG << "FabberIoNewimage::Setting voxel coordinates from extent" << endl;
//int nx = vol.xsize();
//int ny = vol.ysize();
//int nz = vol.zsize();
	volume4D<float> coordvol(nx, ny, nz, 3);
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			for (int k = 0; k < nz; k++)
			{
				ColumnVector vcoord(3);
				vcoord << i << j << k;
				coordvol.setvoxelts(vcoord, i, j, k);
			}
		}
	}

	if (m_have_mask)
	{
		m_coords = coordvol.matrix(m_mask);
	}
	else
	{
		m_coords = coordvol.matrix();
	}
	m_have_coords = true;
}

void FabberIoNewimage::Clear()
{
	m_coords = Matrix();
	m_mask = volume<float>();
	m_have_mask = m_have_coords = false;
}
