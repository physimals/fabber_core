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

enum VoxelDataType
{
  VDT_SCALAR, VDT_MVN
};

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
 * Basic data load/save interface
 */
class FabberIo
{
public:
	virtual ~FabberIo()
	{
	}
	virtual void LoadMask(std::string filename) = 0;
	virtual NEWMAT::Matrix LoadVoxelData(std::string filename) = 0;
	virtual void SaveVoxelData(NEWMAT::Matrix &data, std::vector<int> extent, std::string filename, VoxelDataType data_type) = 0;
	virtual NEWMAT::Matrix GetVoxelCoords() = 0;

};

#endif /* __FABBER_IO_H */
