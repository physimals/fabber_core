/**
 * fabber_io_newimage.h
 *
 * Loading and saving of data to files. This is the
 * only part of Fabber which depends on NEWIMAGE
 * and can be disabled if compiling for an environment
 * which has its own image load/save facilities.
 *
 * Copyright (C) 2007-2008 University of Oxford
 */

#ifndef __FABBER_IO_NEWIMAGE_H
#define __FABBER_IO_NEWIMAGE_H

#ifndef NO_NEWIMAGE

#include "fabber_io.h"

#include "newimage/newimage.h"
#include "newmat.h"

#include <string>
#include <vector>

/**
 * IO module which uses NEWIMAGE to load NIFTII files
 */
class FabberIoNewimage : public FabberIoMemory
{
public:
	FabberIoNewimage();
	void Initialize(FabberRunData &rundata);
	const NEWMAT::Matrix &GetVoxelData(std::string key);
	void SaveVoxelData(NEWMAT::Matrix &data, std::string filename, VoxelDataType data_type);

private:
	void SetVoxelCoordsFromExtent(int nx, int ny, int nz);
	NEWIMAGE::volume<float> m_mask;
	bool m_have_mask;
};

#endif /* NO_NEWIMAGE */

#endif /* __FABBER_IO_NEWIMAGE_H */

