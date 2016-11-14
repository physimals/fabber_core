/**
 * fabber_io.h
 *
 * Interface for loading and saving of data to files.
 * This is not included in the basic Fabber library
 * so it can easily be used in an environment
 * which has its own image load/save facilities.
 *
 * Copyright (C) 2007-2008 University of Oxford
 */

#ifndef __FABBER_IO_H
#define __FABBER_IO_H

#include "newmat.h"

#include <string>
#include <vector>
#include <stdexcept>
#include <map>

/**
 * Type of voxel data.
 *
 * So far we only have scalar (double) data and the MVN structure
 * which associates data to describe a symmetric matrix with each
 * voxel. The size of this data depends on the number of model
 * and noise parameters in a non-trivial way.
 */
enum VoxelDataType
{
	VDT_SCALAR, VDT_MVN
};

// Ugly forward declaration reflacting a circular dependency between
// FabberIo and FabberRunData
class FabberRunData;

/**
 * Thrown when data cannot be loaded from the specified location
 */
class DataLoadError: public runtime_error
{
public:
	DataLoadError(std::string filename) :
			runtime_error("Load error: " + filename)
	{
	}
};

/**
 * Thrown when data cannot be loaded from the specified location
 */
class DataSaveError: public runtime_error
{
public:
	DataSaveError(std::string filename) :
			runtime_error("Save error: " + filename)
	{
	}
};

/**
 * Abstract interface for loading and saving of voxel data
 */
class FabberIo
{
public:
	virtual ~FabberIo()
	{
	}
	/**
	 * Initialize the I/O module.
	 *
	 * Loading/saving modules should be initialized with the run data before use. This enables
	 * them to act on module-specific options (e.g. mask image for the NEWIMAGE
	 * based loader).
	 *
	 * FabberRunData will call this method as part of the Run() method before any data
	 * requests are made.
	 */
	virtual void Initialize(FabberRunData &rundata) = 0;

	/**
	 * Get voxel data for a given key.
	 *
	 * Depending on the implementation, the key may be interpreted as a filename
	 * or as a key to some other storage structure.
	 */
	virtual const NEWMAT::Matrix &GetVoxelData(std::string key) = 0;

	/**
	 * Get the list of voxel co-ordinates
	 *
	 * This method should function regardless of whether or not any actual voxel
	 * data has been loaded yet.
	 */
	virtual const NEWMAT::Matrix &GetVoxelCoords() = 0;

	/**
	 * Save voxel data.
	 *
	 * This will be called for output voxel data that Fabber has been configured to save.
	 * The data key may be interpreted as part of a filename but will not include any
	 * implementation-specific details such as extension or path.
	 */
	virtual void SaveVoxelData(NEWMAT::Matrix &data, std::string filename, VoxelDataType data_type) = 0;

	/**
	 * Get the data extent.
	 *
	 * FIXME dims is not yet implemented
	 *
	 * @param extent will be set to a list of the number of voxels in the x, y, z dimensions
	 * @param dims will be set to a list of the mm physical sizes of each voxel in the x,y, z dimensions,
	 *             if these are available. If not, they will be set equal to extent.
	 */
	virtual void GetExtent(vector<int> &extent, vector<float> &dims) = 0;

	/**
	 * Clear the named voxel data from this IO module.
	 *
	 * If key is not specified, clear all voxel data.
	 *
	 * Internal data set in the Initialize method should not be cleared, only data which
	 * affects LoadVoxelData and SaveVoxelData. GetVoxelCoords should not be affected
	 * by a call to Clear()
	 *
	 * @param key Voxel data to clear, or if not specified, clear all voxel data.
	 */
	virtual void ClearVoxelData(string key="") = 0;
};

/**
 * Load and save data to/from memory
 */
class FabberIoMemory: public FabberIo
{
public:
	FabberIoMemory();
	virtual void Initialize(FabberRunData &rundata);
	virtual const NEWMAT::Matrix &GetVoxelData(std::string key);
	virtual const NEWMAT::Matrix &GetVoxelCoords();
	virtual void SaveVoxelData(NEWMAT::Matrix &data, std::string key, VoxelDataType data_type);
	virtual void GetExtent(vector<int> &extent, vector<float> &dims);
	virtual void ClearVoxelData(string key="");

	/**
	 * Set named voxel data
	 *
	 * @param key Name identifying the voxel data required
	 * @param data an NxM matrix where each column contains the data for a single
	 *        voxel. The rows may contain the time series of data for that voxel
	 *        however they might be used for other purposes, e.g. storing the mean
	 *        of each parameter for that voxel.
	 * @throw If number of columns in data is not equal to the number of voxels
	 */
	virtual void SetVoxelData(string key, const NEWMAT::Matrix &data);

	/**
	 * Set the voxel co-ordinates
	 *
	 * @param coords an Nx3 matrix where each column is a voxel and the rows
	 *         are the xyz co-ords of the voxel. The co-ordinates are
	 *         grid positions (integers), not physical co-ordiantes (mm)
	 */
	void SetVoxelCoords(const NEWMAT::Matrix &coords);


protected:
	void CheckSize(std::string key, const NEWMAT::Matrix &mat);

	int m_nvoxels;
	NEWMAT::Matrix m_voxel_coords;
	bool m_have_coords;
	std::map<std::string, NEWMAT::Matrix> m_voxel_data;
	std::vector<int> m_extent;
	std::vector<float> m_dims;
};

#endif /* __FABBER_IO_H */
