/**
 * fabber_rundata_recxml.h
 *
 * Extends FabberRunData with loading and saving of data to
 * RECXML files. This is a proprietory Phillips format 
 * consisting of a header file in XML format and a binary
 * data file
 *
 * Martin Craig
 *
 * Copyright (C) 2019 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "rundata_newimage.h"

#include <newmat.h>

#include <string>

struct RecxmlImage
{   
    std::map<std::string, std::string> global_attrs;
    std::vector<std::map<std::string, std::string> > vol_attrs;
};

/**
 * Run data which loads/saves to/from RECXML files
 */
class FabberRunDataRecxml : public FabberRunDataNewimage
{
public:
    FabberRunDataRecxml(bool compat_options = true);

    void SetExtentFromData();
    const NEWMAT::Matrix &LoadVoxelData(const std::string &filename);
    virtual void SaveVoxelData(
        const std::string &filename, NEWMAT::Matrix &data, VoxelDataType data_type = VDT_SCALAR);

    // For metadata parser
    bool m_parser_inimage;
    std::string m_parser_currdata;
    std::string m_parser_attrname;
    std::string m_parser_chars;
    std::map<std::string, RecxmlImage> m_recxml_data;

private:
    void LoadMetadata(const std::string &filename);
    void LoadBinaryData(const std::string &filename);
};
