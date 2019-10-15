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

FabberRunDataRecxml::FabberRunDataRecxml(bool compat_options)
    : FabberRunDataNewimage(compat_options)
{
}

void FabberRunDataRecxml::SetExtentFromData()
{
    // FIXME hardcoded for now
    SetCoordsFromExtent(96, 96, 11);
}

const Matrix &FabberRunDataRecxml::LoadVoxelData(const std::string &filename)
{
    if (m_voxel_data.find(filename) == m_voxel_data.end())
    {
        LOG << "FabberRunDataRecxml::Loading data from '" + filename << "'" << endl;

        volume4D<float> vol;
        m_parser_currdata = filename;
        try
        {
            // Load the metadata file
            LoadMetadata(filename + ".xml");

            // Load the binary data into a NEWIMAGE volume
            LoadBinaryData(filename + ".rec");
            if (!m_have_mask)
            {
                // We need a mask volume so that when we save we can make sure
                // the image properties are set consistently with the source data
                m_mask = vol[0];
                m_mask = 1;
                m_have_mask = true;
            }
        }
        catch (...)
        {
            throw DataNotFound(filename, "Error loading file");
        }
        //DumpVolumeInfo4D(vol, LOG);

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
    }

    return m_voxel_data[filename];
}

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
            cerr << "Image ";
            // Currently inside an image record - add attribute to the last image in the sequence
            int last_vol = recxml.m_recxml_data[recxml.m_parser_currdata].vol_attrs.size() - 1;
            recxml.m_recxml_data[recxml.m_parser_currdata].vol_attrs[last_vol][recxml.m_parser_attrname] = recxml.m_parser_chars;
        }
        else 
        {
            cerr << "Global ";
            // Not currently inside an image record - add attribute to global list
            recxml.m_recxml_data[recxml.m_parser_currdata].global_attrs[recxml.m_parser_attrname] = recxml.m_parser_chars;
        }
        cerr << "Attribute: " << recxml.m_parser_attrname << "=" << recxml.m_parser_chars << endl;
    }
    recxml.m_parser_chars = "";
}

static void characters_handler(void * ctx, const xmlChar * ch, int len)
{
    FabberRunDataRecxml &recxml = *( static_cast<FabberRunDataRecxml *>( ctx ) );
    std::string chars((const char *)ch, (const char *)ch+len);
    recxml.m_parser_chars += chars;
}

void FabberRunDataRecxml::LoadMetadata(const std::string &filename)
{
   LIBXML_TEST_VERSION

   xmlSAXHandler saxHandler; // See http://xmlsoft.org/html/libxml-tree.html#xmlSAXHandler
   memset( &saxHandler, 0, sizeof(saxHandler) );
   // Using xmlSAXVersion( &saxHandler, 2 ) generate crash as it sets plenty of other pointers...
   saxHandler.initialized = XML_SAX2_MAGIC;  // so we do this to force parsing as SAX2.
   saxHandler.startElementNs = &startElementNs;
   saxHandler.endElementNs = &endElementNs;
   saxHandler.characters = &characters_handler;

   m_parser_inimage = false;
   int result = xmlSAXUserParseFile(&saxHandler, this, filename.c_str());
   xmlCleanupParser(); 
   if (result != 0)
   {
      throw std::runtime_error("Failed to parse document");
   }
}

void FabberRunDataRecxml::LoadBinaryData(const std::string &filename)
{
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
