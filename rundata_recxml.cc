/**
 * rundata_recxml.cc
 *
 * Martin Craig
 *
 * Copyright (C) 2019 University of Oxford  */

/*  CCOPYRIGHT */

#include "rundata_recxml.h"

#include "easylog.h"
#include "rundata.h"

#include <newimage/newimage.h>
#include <newimage/newimageio.h>
#include <newmat.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

#include <ostream>
#include <string>
#include <vector>

using namespace std;
using NEWMAT::Matrix;
using NEWIMAGE::volume;

#define MAX_BYTES_PER_VALUE 8

/**
 * SAX callback for the start of an XML element
 */
static void startElementNs( void * ctx, 
                            const xmlChar * localname, 
                            const xmlChar * prefix, 
                            const xmlChar * URI, 
                            int nb_namespaces, 
                            const xmlChar ** namespaces, 
                            int nb_attributes, 
                            int nb_defaulted, 
                            const xmlChar ** attributes )
{
    FabberRunDataRecxml &recxml = *( static_cast<FabberRunDataRecxml *>( ctx ) );
    if (string((const char *)localname) == "Image_Info") 
    {
        // Start of an image record - add an attribute map for it
        recxml.m_parser_inimage = true;
        recxml.m_recxml_data[recxml.m_parser_currdata].vol_attrs.push_back(std::map<std::string, std::string>());
    }
    else if (string((const char *)localname) == "Attribute")
    {
        unsigned int index = 0;
        for ( int indexAttribute = 0; 
            indexAttribute < nb_attributes; 
            ++indexAttribute, index += 5 )
        {
            const xmlChar *localname = attributes[index];
            const xmlChar *valueBegin = attributes[index+3];
            const xmlChar *valueEnd = attributes[index+4];
            std::string value( (const char *)valueBegin, (const char *)valueEnd );
            if (string((const char *)localname) == "Name")
            {
                recxml.m_parser_attrname = value;
            }
        }
    }
}

/**
 * SAX callback for the end of an XML element
 */
static void endElementNs( void * ctx, 
                          const xmlChar * localname, 
                          const xmlChar * prefix, 
                          const xmlChar * URI )
{
    FabberRunDataRecxml &recxml = *( static_cast<FabberRunDataRecxml *>( ctx ) );
    if (string((const char *)localname) == "Image_Info") 
    {
        recxml.m_parser_inimage = false;
    }
    else if (string((const char *)localname) == "Attribute")
    {
        if (recxml.m_parser_inimage) 
        {
            //cerr << "Image ";
            // Currently inside an image record - add attribute to the last image in the sequence
            int last_vol = recxml.m_recxml_data[recxml.m_parser_currdata].vol_attrs.size() - 1;
            recxml.m_recxml_data[recxml.m_parser_currdata].vol_attrs[last_vol][recxml.m_parser_attrname] = recxml.m_parser_chars;
        }
        else 
        {
            //cerr << "Global ";
            // Not currently inside an image record - add attribute to global list
            recxml.m_recxml_data[recxml.m_parser_currdata].global_attrs[recxml.m_parser_attrname] = recxml.m_parser_chars;
        }
        //cerr << "Attribute: " << recxml.m_parser_attrname << "=" << recxml.m_parser_chars << endl;
    }
    recxml.m_parser_chars = "";
}

/**
 * SAX callback for between-element characters
 */
static void characters_handler(void * ctx, const xmlChar * ch, int len)
{
    FabberRunDataRecxml &recxml = *( static_cast<FabberRunDataRecxml *>( ctx ) );
    std::string chars((const char *)ch, (const char *)ch+len);
    recxml.m_parser_chars += chars;
}

FabberRunDataRecxml::FabberRunDataRecxml(bool compat_options)
    : FabberRunDataNewimage(compat_options)
{
}

void FabberRunDataRecxml::SetExtentFromData()
{
    string mask_fname = GetStringDefault("mask", "");
    m_have_mask = (mask_fname != "");

    if (m_have_mask)
    {
        LOG << "FabberRunDataNewimage::Loading mask data from '" + mask_fname << "'" << endl;
        LoadRecxml(m_mask, mask_fname);
        m_mask.binarise(1e-16, m_mask.max() + 1, NEWIMAGE::exclusive);
        DumpVolumeInfo(m_mask);
        SetCoordsFromExtent(m_mask.xsize(), m_mask.ysize(), m_mask.zsize());
    }
    else
    {
        // Make sure the coords are loaded from the main data even if we don't
        // have a mask, and that the reference volume is initialized
        LOG << "FabberRunDataNewimage::No mask, using data for extent" << endl;
        string data_fname = GetStringDefault("data", GetStringDefault("data1", ""));
        volume4D<float> main_vol;
        LoadRecxml(main_vol, data_fname);
        SetCoordsFromExtent(main_vol.xsize(), main_vol.ysize(), main_vol.zsize());
        m_voxel_data[data_fname] = main_vol.matrix();
        m_mask = main_vol[0];
        m_mask = 1;
        m_have_mask = true;
    }
}

const Matrix &FabberRunDataRecxml::LoadVoxelData(const std::string &filename)
{
    if (m_voxel_data.find(filename) == m_voxel_data.end())
    {
        LOG << "FabberRunDataRecxml::Loading data from '" + filename << "'" << endl;

        volume4D<float> vol;
        try
        {
            LoadRecxml(vol, filename);
        }
        catch (...)
        {
            throw DataNotFound(filename, "Error loading file");
        }
        DumpVolumeInfo4D(vol);

        try
        {
            if (m_have_mask)
            {
                LOG << "FabberRunDataRecxml::Applying mask to data..." << endl;
                m_voxel_data[filename] = vol.matrix(m_mask);
            }
            else
            {
                m_voxel_data[filename] = vol.matrix();
            }
        }
        catch (exception &e)
        {
            LOG << "NEWMAT error while applying mask... Most likely a dimension mismatch. ***\n";
            throw;
        }

        SaveVoxelData(filename + "_temp", m_voxel_data[filename]);
    }

    return m_voxel_data[filename];
}

void FabberRunDataRecxml::LoadRecxml(volume4D<float> &vol, const std::string &filename)
{
    // Parse XML header using libxml SAX interface
    // FIXME should have a parser class and keep the messy flags out of the rundata class
    LIBXML_TEST_VERSION
    xmlSAXHandler saxHandler;
    memset( &saxHandler, 0, sizeof(saxHandler) );
    saxHandler.initialized = XML_SAX2_MAGIC;
    saxHandler.startElementNs = &startElementNs;
    saxHandler.endElementNs = &endElementNs;
    saxHandler.characters = &characters_handler;

    m_parser_inimage = false;
    m_parser_currdata = filename;
    int result = xmlSAXUserParseFile(&saxHandler, this, (filename + ".xml").c_str());
    xmlCleanupParser(); 
    if (result != 0)
    {
       throw DataNotFound(filename, "Failed to parse XML header");
    }

    // Extract the XY resolution from the first image - FIXME check all the same
    int rx = convertTo<int>(m_recxml_data[m_parser_currdata].vol_attrs[0]["Resolution X"]);
    int ry = convertTo<int>(m_recxml_data[m_parser_currdata].vol_attrs[0]["Resolution Y"]);
    int px = convertTo<int>(m_recxml_data[m_parser_currdata].vol_attrs[0]["Pixel Size"]);
    LOG << "FabberRunDataRecxml::Resolution: " << rx << "x" << ry << " at " << px << " bits/pixel" << endl;
    int bytes_per_value = px / 8;
    if (bytes_per_value > MAX_BYTES_PER_VALUE)
    {
        throw DataNotFound(filename, "Data bits / pixel too large: " + stringify(px));
    }

    // Find the max number of slices and phases - this gives the Z and T resolution
    int n_slices = 0;
    int n_phases = 0;
    for (size_t idx=0; idx<m_recxml_data[m_parser_currdata].vol_attrs.size(); idx++) 
    {
        int slice = convertTo<int>(m_recxml_data[m_parser_currdata].vol_attrs[idx]["Slice"]);
        int phase = convertTo<int>(m_recxml_data[m_parser_currdata].vol_attrs[idx]["Phase"]);
        if (slice > n_slices) n_slices = slice;
        if (phase > n_phases) n_phases = phase;
    }
    vol.reinitialize(rx, ry, n_slices, n_phases);

    // Now read the binary file and put the data in the appropriate place
    FILE *fp = fopen((filename + ".rec").c_str(), "rb");
    if (!fp)
    {
        throw DataNotFound(filename, "Failed to open binary image data file");
    }

    unsigned char bytes[MAX_BYTES_PER_VALUE];
    for (size_t idx=0; idx<m_recxml_data[m_parser_currdata].vol_attrs.size(); idx++) 
    {
        int slice = convertTo<int>(m_recxml_data[m_parser_currdata].vol_attrs[idx]["Slice"]);
        int phase = convertTo<int>(m_recxml_data[m_parser_currdata].vol_attrs[idx]["Phase"]);
        for (int x=0; x<rx; x++) 
        {
            for (int y=0; y<ry; y++) 
            {
                // Reading values one at a time - possibly rather slow but not causing problems
                // so far
                int values_read = fread(bytes, bytes_per_value, 1, fp);
                if (values_read != 1)
                {
                    fclose(fp);
                    throw DataNotFound(filename, "Unexpected end-of-file reading binary image data");
                }

                // REC files are apparently always little-endian - this is to make sure we 
                // read them as such regardless of platform endianness FIXME can we be sure
                // data is always integers?
                unsigned int value = 0;
                for (int byte=0; byte<bytes_per_value; byte++) 
                {
                    value |= bytes[byte] << (8*byte);
                }
                vol.value(x, y, slice-1, phase-1) = float(value);
            }
        }
    }
    fclose(fp);
}

void FabberRunDataRecxml::SaveVoxelData(
    const std::string &filename, NEWMAT::Matrix &data, VoxelDataType data_type)
{
    LOG << "FabberRunDataRecxml::Saving to nifti: " << filename << endl;
    int nifti_intent_code;
    switch (data_type)
    {
    case VDT_MVN:
        nifti_intent_code = NIFTI_INTENT_SYMMATRIX;
        break;
    default:
        nifti_intent_code = NIFTI_INTENT_NONE;
    }

    int data_size = data.Nrows();
    volume4D<float> output(m_extent[0], m_extent[1], m_extent[2], data_size);
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

    if (filename[0] == '/')
    {
        // Absolute path
        save_volume4D(output, filename);
    }
    else
    {
        // Relative path
        string filepath = GetOutputDir() + "/" + filename;
        save_volume4D(output, filepath);
    }
}
