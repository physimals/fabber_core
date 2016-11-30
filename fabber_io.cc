#include "fabber_io.h"

#include "dataset.h"
#include "easylog.h"

#include "newmat.h"

#include <string>
#include <map>

using namespace std;
using NEWMAT::Matrix;

FabberIoMemory::FabberIoMemory() :
		m_nvoxels(-1), m_have_coords(false)
{
}

void FabberIoMemory::Initialize(FabberRunData &rundata)
{
}

const NEWMAT::Matrix &FabberIoMemory::GetVoxelData(std::string key)
{
	if (m_voxel_data.count(key) == 0)
	{
		throw DataNotFound(key);
	}
	else
	{
		return m_voxel_data.find(key)->second;
	}
}

void FabberIoMemory::SaveVoxelData(NEWMAT::Matrix &data, std::string filename, VoxelDataType data_type)
{
	LOG << "FabberIoMemory::Saving to memory: " << filename << endl;
	// FIXME what should we do with data_type?
	SetVoxelData(filename, data);
}

const NEWMAT::Matrix &FabberIoMemory::GetVoxelCoords()
{
	// Throw exception if we have no co-ordinates - i.e. user has not
	// called SetVoxelCoords yet.
	if (m_have_coords)
	{
		return m_voxel_coords;
	}
	else
	{
		throw DataNotFound("voxel coordinates");
	}
}

void FabberIoMemory::ClearVoxelData(string key)
{
	if (key != "")
	{
		m_voxel_data.erase(key);
	}
	else
	{
		m_voxel_data.clear();

	}

	if ((m_voxel_data.size() == 0) && !m_have_coords)
	{
		m_nvoxels = -1;
	}
}

void FabberIoMemory::SetVoxelData(string key, const NEWMAT::Matrix &data)
{
	CheckSize(key, data);
	m_voxel_data[key] = data;
}

void FabberIoMemory::SetVoxelCoords(const NEWMAT::Matrix &coords)
{
	if (m_voxel_data.size() == 0) {
		// No data, so make sure coords overrides existing data
		m_nvoxels = -1;
		m_have_coords = false;
	}
	// This will set m_nvoxels if we don't already have data
	CheckSize("coords", coords);

	// We assume 3D coordinates. Fabber could work for different
	// numbers of dimensions but would require extensive refactoring
	if (coords.Nrows() != 3)
	{
		throw Invalid_option("Co-ordinates must be 3 dimensional");
	}

	// FIXME we assume coords will not be negative
	m_extent.resize(3);
	m_extent[0] = coords.Row(1).Maximum() - coords.Row(1).Minimum() + 1;
	m_extent[1] = coords.Row(2).Maximum() - coords.Row(2).Minimum() + 1;
	m_extent[2] = coords.Row(3).Maximum() - coords.Row(3).Minimum() + 1;

	m_dims.resize(3);
	m_dims[0] = 1.0;
	m_dims[1] = 1.0;
	m_dims[2] = 1.0;

	m_voxel_coords = coords;
	m_have_coords = true;
}

void FabberIoMemory::GetExtent(vector<int> &extent, vector<float> &dims)
{
	if (m_have_coords)
	{
		extent = m_extent;
		dims = m_dims;
	}
	else
	{
		throw DataNotFound("voxel coordinates");
	}
}

void FabberIoMemory::CheckSize(std::string key, const NEWMAT::Matrix &mat)
{
	if (m_nvoxels == -1)
	{
		// If this is the first data set provided, it gets to
		// decide the number of voxels
		m_nvoxels = mat.Ncols();
	}
	else
	{
		// Otherwise, check data provided has same number of voxels
		// as previous data
		if (mat.Ncols() != m_nvoxels)
		{
			throw Invalid_option(
					"Per-voxel matrix " + key + " is incorrect size (cols=" + stringify(mat.Ncols()) + " should be "
							+ stringify(m_nvoxels) + ")");
		}
	}
}
