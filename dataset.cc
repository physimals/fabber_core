/*  dataset.cc - Data-loading class for FABBER

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "easylog.h"
#include "dataset.h"

using namespace MISCMATHS;
using namespace std;

#ifndef __FABBER_LIBRARYONLY 
using namespace NEWIMAGE;

void DumpVolumeInfo(const volume4D<float>& info, string indent = "", ostream& out = LOG)
{
  Tracer_Plus tr("DumpVolumeInfo");
  LOG << indent << "Dimensions: x=" << info.xsize() << ", y=" << info.ysize() 
      << ", z=" << info.zsize() << ", vols=" << info.tsize() << endl;
  LOG << indent << "Voxel size: x=" << info.xdim() << "mm, y=" << info.ydim() 
      << "mm, z=" << info.zdim() << "mm, TR=" << info.tdim() << " sec\n";
  LOG << indent << "Intents: " << info.intent_code() << ", " << info.intent_param(1)
      << ", " << info.intent_param(2) << ", " << info.intent_param(3) << endl;
}

void DumpVolumeInfo(const volume<float>& info, string indent = "", ostream& out = LOG)
{
  Tracer_Plus tr("DumpVolumeInfo");
  LOG << indent << "Dimensions: x=" << info.xsize() << ", y=" << info.ysize() 
      << ", z=" << info.zsize() << ", vols=1" << endl;
  LOG << indent << "Voxel size: x=" << info.xdim() << "mm, y=" << info.ydim() 
      << "mm, z=" << info.zdim() << "mm, TR=1" << " sec\n";
  LOG << indent << "Intents: " << info.intent_code() << ", " << info.intent_param(1)
      << ", " << info.intent_param(2) << ", " << info.intent_param(3) << endl;
}

#endif //__FABBER_LIBRARYONLY

// Inputs: reads various options from args, and loads the input data
// Outputs: masks is set (except with UsingMatrixIO) and voxelData is populated
void DataSet::LoadData(ArgsType& args)
{
  Tracer_Plus tr("LoadData");

  if (EasyOptions::UsingMatrixIO())
    {
      string dataFile = args.Read("data");
      voxelData = EasyOptions::InMatrix(dataFile);

      string voxelCoordsFile = args.ReadWithDefault("voxelCoords","");
      if (voxelCoordsFile != "")
      {
	voxelCoords = EasyOptions::InMatrix(voxelCoordsFile);
	// leave mask undefined
      }
      else
      {
        // voxelCoords is left undefined -- hope it's not used!!
      }

      string suppdataFile = args.ReadWithDefault("suppdata","none");
      if (suppdataFile != "none") {
	voxelSuppData = EasyOptions::InMatrix(suppdataFile);
      }

      return;
    }

#ifdef __FABBER_LIBRARYONLY
  throw Logic_error("NEWIMAGE support was not compiled in, but UsingMatrixIO is false!");
#else
  string dataOrder = args.ReadWithDefault("data-order","interleave");

  // mask
      string maskFile = args.Read("mask");
      LOG_ERR("    Loading mask data from '" + maskFile << "'" << endl);
      read_volume(mask,maskFile);
      DumpVolumeInfo(mask);

      // create coordinates matrix using mask
      int nx = mask.xsize();
      int ny = mask.ysize();
      int nz = mask.zsize();
      volume4D<float> coordvol(nx,ny,nz,3);
      for (int i=0; i<nx; i++) { for (int j=0; j<ny; j++) { for (int k=0; k<nz; k++) {
	    ColumnVector vcoord(3);
	    vcoord << i << j << k;
	    coordvol.setvoxelts(vcoord,i,j,k);
	  } } }
      voxelCoords=coordvol.matrix(mask);
      mask.binarise(1e-16,mask.max()+1,exclusive);

      // supplementary data
      string suppdataFile = args.ReadWithDefault("suppdata","none");
      if (suppdataFile != "none") {
	LOG_ERR("    Loading supplementary data from '" + suppdataFile << "'" << endl);
	volume4D<float> suppdata;
	read_volume4D(suppdata,suppdataFile);
	DumpVolumeInfo(suppdata);

	LOG << "     Applying mask to supplementary data..." << endl;
	try {
	  voxelSuppData=suppdata.matrix(mask);
	} catch (exception) {
	  LOG_ERR("\n*** NEWMAT error while thresholding SUPPLEMENTARY time-series... "
		  << "Most likely a dimension mismatch. ***\n");
        DumpVolumeInfo(suppdata);
        LOG_ERR("Mask:\n");
        DumpVolumeInfo(mask);
        LOG_ERR("\nThis is fatal... rethrowing exception.\n");
        throw;
      }

      }

  if (dataOrder == "singlefile")
    {
      LOG << "  Loading all data from a single file..." << endl;
      string dataFile = args.Read("data");

      // Load the files.  Note: these functions don't throw errors if file doesn't exist --
      // they just crash.  Hence the detailed logging before we try anything.

      LOG_ERR("    Loading data from '" << dataFile << "'" << endl);
      volume4D<float> data;
      read_volume4D(data,dataFile);
      DumpVolumeInfo(data);

      

      LOG << "    Applying mask to data..." << endl;
      //mask.binarise(1e-16,mask.max()+1,exclusive);
      // threshold using mask:
      try {
	voxelData=data.matrix(mask);
      } catch (Exception) {
        LOG_ERR("\n*** NEWMAT error while thresholding time-series... "
                << "Most likely a dimension mismatch. ***\n");
        DumpVolumeInfo(data);
        LOG_ERR("Mask:\n");
        DumpVolumeInfo(mask);
        LOG_ERR("\nThis is fatal... rethrowing exception.\n");
        throw;
      }

    }
  else if (dataOrder == "interleave" || dataOrder == "concatenate")
    {
      LOG << "  Loading data from multiple files..." << endl;
      
      vector<volume4D<float> > dataSets;
      vector<Matrix> dataSetsM;
      int nTimes = -1;
      while (true)
	{
	  int N = dataSets.size() + 1;
	  string datafile = args.ReadWithDefault("data"+stringify(N), "stop!");
	  if (datafile == "stop!") break;

	  // Load the files.  Note: these functions don't throw errors if file doesn't exist --
	  // they just crash.  Hence the detailed logging before we try anything.
	  LOG_ERR("    Loading " << "data"+stringify(N) << " from '" << datafile << "'" << endl);
          volume4D<float> temp;
	  read_volume4D(temp,datafile);
          dataSets.push_back(temp);
	  DumpVolumeInfo(dataSets.back(), "      ");

	  if (nTimes == -1) 
	    nTimes = dataSets[0].tsize();
	  else if ( ( nTimes != dataSets.back().tsize() ) & dataOrder == "interleave")
	    // data sets only strictly need same number of time points if they are to be interleaved
	    throw Invalid_option("Data sets must all have the same number of time points");
	}

      int nSets = dataSets.size();
      if (nSets < 1)
	throw Invalid_option("At least one data file is required: --data1=<file1> [--data2=<file2> [...]]\n");      

      //string maskFile = args.Read("mask");
 
      //LOG_ERR("    Loading mask data from '" + maskFile << "'" << endl);
      //read_volume(mask,maskFile);
      //DumpVolumeInfo(mask);

      LOG << "    Applying mask to all data sets..." << endl;
      //mask.binarise(1e-16,mask.max()+1,exclusive);
      // threshold using mask:
      for (int i = 0; i < nSets; i++)
	{
	  try {
	    dataSetsM.push_back(dataSets[i].matrix(mask));

            // Note that the above doesn't catch all dimension mismatches..
            // If the mask is smaller (in the z-dir, at least) than the data,
            // it doesn't seem to raise any exception.

            if (dataSets[i].xsize() != mask.xsize())
              LOG_ERR("Warning: nonfatal dimension mismatch in x!\n");
            if (dataSets[i].ysize() != mask.ysize())
              LOG_ERR("Warning: nonfatal dimension mismatch in y!\n");
            if (dataSets[i].zsize() != mask.zsize())
              LOG_ERR("Warning: nonfatal dimension mismatch in z!\n");

	  } catch (Exception) {
	    LOG_ERR("\n*** NEWMAT error while thresholding time-series... "
		    << "Most likely a dimension mismatch (more details in logfile) ***\n");
	    LOG << "Data set " << i+1 << ":\n"; 
	    DumpVolumeInfo(dataSets.at(i));
	    LOG << "Mask:\n"; 
	    DumpVolumeInfo(mask);
	    LOG_ERR("\nThis is fatal... rethrowing exception.\n");
	    throw;
	  }
	}
      
      if (dataOrder == "interleave")
        {
          LOG << "    Combining data into one big matrix by interleaving..." << endl;
          // Interleave:
          voxelData.ReSize(nTimes * nSets, dataSetsM[0].Ncols());
          for (int i = 0; i < nTimes; i++)
            {
              for (int j = 0; j < nSets; j++)
	        {
	          voxelData.Row(nSets*i+j+1) = dataSetsM.at(j).Row(i+1);
	        }
	    }
        }
      else
        {
          LOG << "    Combining data into one big matrix by concatenating..." << endl;
          // Concatenate:
          voxelData = dataSetsM.at(0);
          for (unsigned j = 1; j < dataSetsM.size(); j++)
            voxelData &= dataSetsM.at(j);
        }
      
      LOG << "    Done loading data, size = " 
	  << voxelData.Nrows() << " timepoints by "
	  << voxelData.Ncols() << " voxels" << endl;
    }
  else
    throw Invalid_option(("Unrecognized --dataorder: " + dataOrder + " (try interleave or singlefile)").c_str());
#endif //__FABBER_LIBRARYONLY
}

