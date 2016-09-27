/*  dataset.cc - Data-loading class for FABBER

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "newmat.h"

#include "utils/tracer_plus.h"
#include "easylog.h"
#include "setup.h"
#include "fwdmodel.h"
#include "inference.h"
#include "dataset.h"

using namespace std;

void PercentProgressCheck::operator()(int voxel, int nVoxels)
{
	int percent = (100 * voxel) / nVoxels;
	if (percent / 10 > m_last)
	{
		cout << "\b\b\b";
		m_last = percent / 10;
		if (m_last == 0)
			cout << " ";
		cout << m_last * 10 << "%" << flush;
		if (m_last == 10)
			cout << endl;
	}
}

FabberRunData::FabberRunData() :
		m_save_files(false), m_progress(0), m_have_coords(false), m_nvoxels(-1)
{
	FabberSetup::SetupDefaults();
}

FabberRunData::~FabberRunData()
{
	FabberSetup::Destroy();
}

void FabberRunData::LogParams()
{
	map<string, string>::iterator iter;
	for (iter = m_params.begin(); iter != m_params.end(); iter++)
	{
		LOG << "FabberRunData::Parameter " << iter->first << "=" << iter->second << endl;
	}
}
void FabberRunData::Run()
{
	Tracer_Plus tr2("FabberRunData::Run");

	LogParams();

	//Set the forward model
	std::auto_ptr<FwdModel> fwd_model(FwdModel::NewFromName(GetString("model")));
	fwd_model->Initialize(*this);
	assert(fwd_model->NumParams() > 0);
	LOG << "FabberRunData::Forward Model version " << fwd_model->ModelVersion() << endl;

	//Set the inference technique (and pass in the model)
	std::auto_ptr<InferenceTechnique> infer(InferenceTechnique::NewFromName(GetString("method")));
	infer->Initialize(fwd_model.get(), *this);

	// Arguments should all have been used by now, so complain if there's anything left.
	// FIXME ineffective at present
	CheckEmpty();

	// Calculations
	Progress(0, m_nvoxels);
	infer->DoCalculations(*this);
	Progress(m_nvoxels, m_nvoxels);
	infer->SaveResults(*this);
}

static string trim(string const& str)
{
    if(str.empty())
        return str;

    size_t firstScan = str.find_first_not_of(' ');
    size_t first     = firstScan == string::npos ? str.length() : firstScan;
    size_t last      = str.find_last_not_of(' ');
    return str.substr(first, last-first+1);
}

void FabberRunData::ParseParamFile(const string filename)
{
	Tracer_Plus tr("FabberRunData::ParseParamFile");
	ifstream is(filename.c_str());
	if (!is.good())
	{
		LOG << "FabberRunData::Couldn't read input options file: '" << filename << "'";
		throw Invalid_option("Couldn't read input options file:" + filename);
	}
	char c;
	string param;
	while (is.good())
	{
		string input;
		std::getline(is, input);
		input = trim(input);
		if (input.size() > 0)
		{
			if (input[0] == '#')
				continue;
			else
			{
				AddKeyEqualsValue(input, true);
			}
		}
	}
}

void FabberRunData::ParseOldStyleParamFile(const string filename)
{
	Tracer_Plus tr("FabberRunData::ParseOldStyleParamFile");
	ifstream is(filename.c_str());
	if (!is.good())
	{
		LOG << "FabberRunData::Couldn't read -@ input file: '" << filename << "'";
		throw Invalid_option("Couldn't read input file: -@ " + filename);
	}
	char c;
	string param;
	while (is.good())
	{
		is.get(c);
		if (!isspace(c))
			param += c;
		else if (param == "")
		{
		} // repeated whitespace, so do nothing
		else if (string(param, 0, 2) == "--") // we have an option
		{
			AddKeyEqualsValue(string(param, 2, string::npos));
			param = "";
		}
		else if (string(param, 0, 1) == "#") // comment
		{
			// discard this word and the rest of the line.
			param = "";
			while (is.good() && c != '\n')
				is.get(c);
		}
		else if (string(param, 0, 2) == "-@")
		{
			throw Invalid_option("Sorry, at the moment you can only use -@ on the command line (no recursion).\n");
		}
		else
		{
			throw Invalid_option("Invalid option '" + param + "' found in file '" + filename + "'\n");
		}
	}
}

void FabberRunData::Parse(int argc, char** argv)
{
	Tracer_Plus tr("FabberRunData::Parse");
	m_save_files = true; // FIXME hack for CL tool

	m_params[""] = argv[0];
	for (int a = 1; a < argc; a++)
	{
		if (string(argv[a]) == "-f")
		{
			// FIXME what if no file specified?
			ParseParamFile(argv[++a]);
			++a;
		}
		else if (string(argv[a], 0, 2) == "--")
		{
			string key = argv[a] + 2; // skip the "--"
			AddKeyEqualsValue(key);
		}
		else if (string(argv[a]) == "-@")
		{
			ParseOldStyleParamFile(argv[++a]);
		}
		else
		{
			LOG << "FabberRunData::Option doesn't begin with --: " + string(argv[a]);
			throw Invalid_option("An option doesn't begin with --\n");
		}
	}
}

void FabberRunData::AddKeyEqualsValue(const string exp, bool trim_comments)
{
	string::size_type eqPos = exp.find("=");
	string key = trim(string(exp, 0, eqPos));
	if (m_params.count(key) > 0)
		throw Invalid_option("Duplicated option: '" + key + "'\n");
	else if (eqPos != (exp.npos)) {
		string::size_type end = exp.npos;
		if (trim_comments) end = exp.find("#");
		string value = trim(exp.substr(eqPos + 1, end-(eqPos+1)));
		m_params[key] = value;
	}
	else
	{
		m_params[exp] = "";
	}
}

void FabberRunData::Set(const string key, const string value)
{
	m_params[key] = value;
}

void FabberRunData::Set(const string key, double value)
{
	m_params[key] = stringify(value);
}

void FabberRunData::SetBool(const string key, bool value)
{
	if (value)
		m_params[key] = "";
	else
		m_params.erase(key);
}

string FabberRunData::GetString(const string key)
{
	return GetString(key, "Missing mandatory option: --" + key + "\n");
}

string FabberRunData::GetString(const string key, const string msg)
{
	if (m_params.count(key) == 0)
		throw Invalid_option(msg);

	if (m_params[key] == "")
		throw Invalid_option("No value given for mandatory option: --" + key + "=???");

	// okay, option is valid.  Now remove it.
	string ret = m_params[key];
	//    m_params.erase(key);
	return ret;
}

bool FabberRunData::GetBool(const string key)
{
	if (m_params.count(key) == 0)
		return false;

	if (m_params[key] == "")
	{
		//       m_params.erase(key);
		return true;
	}

	throw Invalid_option("Value should not be given for boolean option --" + key);
}

string FabberRunData::GetStringDefault(const string key, const string def)
{
	if (m_params.count(key) == 0)
		return def;
	if (m_params[key] == "")
		throw Invalid_option("Option requires a value: --" + key + " (or omit, equivalent to --" + key + "=" + def);
	string ret = m_params[key];
	//   m_params.erase(key);
	return ret;
}

void FabberRunData::CheckEmpty()
{
	m_params.erase(""); // not worth complaining about this

	if (m_params.empty())
		return;

	//    string msg = "\nUnused arguments:\n" + stringify(*this);
	//    throw Invalid_option(msg);
}

ostream& operator<<(ostream& out, const FabberRunData& opts)
{
	for (map<string, string>::const_iterator i = opts.m_params.begin(); i != opts.m_params.end(); i++)
	{
		if (i->second == "")
			out << "--" << i->first << endl;
		else
			out << "--" << i->first << "='" << i->second << "'" << endl;
	}
	return out;
}

const NEWMAT::Matrix& FabberRunData::GetMainVoxelData()
{
	// Main voxel data is a bit special because it can
	// come from multiple files
	try
	{
		return GetVoxelData("data");
	} catch (DataNotFound &e)
	{
		GetMainVoxelDataMultiple();
		return GetVoxelData("data");
	}
}

void FabberRunData::SetMainVoxelData(NEWMAT::Matrix &data)
{
	FabberRunData::SetVoxelData("data", data);
}

const NEWMAT::Matrix& FabberRunData::GetVoxelSuppData()
{
	// FIXME Currently Fabber models assume that suppdata will return an empty matrix if none present.
	try
	{
		return GetVoxelData("suppdata");
	} catch (DataNotFound &e)
	{
		return m_empty;
	}
}

void FabberRunData::SetVoxelSuppData(NEWMAT::Matrix &data)
{
	FabberRunData::SetVoxelData("suppdata", data);
}

/**
 * Check matrix is per-voxel, i.e. the number of columns
 * equals the number of voxels
 */
void FabberRunData::CheckSize(std::string key, NEWMAT::Matrix &mat)
{
	// If this is the first data set provided, it gets to
	// decide the number of voxels
	if (m_nvoxels == -1)
	{
		m_nvoxels = mat.Ncols();
	}
	else
	{
		if (mat.Ncols() != m_nvoxels)
		{
			throw Invalid_option(
					"Per-voxel matrix " + key + " is incorrect size (cols=" + stringify(mat.Ncols()) + " should be "
							+ stringify(m_nvoxels) + ")");
		}
	}
}

void FabberRunData::SetVoxelCoords(NEWMAT::Matrix &coords)
{
	// This will also set m_nvoxels if it's the first data set
	CheckSize("coords", coords);
	m_voxelCoords = coords;

	// FIXME we assume 3D coordinates
	if (coords.Nrows() != 3)
	{
		throw Invalid_option("Co-ordinates must be 3 dimensional");
	}

	// FIXME we assume coords will start at 0 in any
	// volume that is output but will not be negative
	m_size.resize(3);
	m_size[0] = coords.Row(1).Maximum() + 1;
	m_size[1] = coords.Row(2).Maximum() + 1;
	m_size[2] = coords.Row(3).Maximum() + 1;

	m_dims.resize(3);
	m_dims[0] = 1.0;
	m_dims[1] = 1.0;
	m_dims[2] = 1.0;

	m_have_coords = true;
}

int FabberRunData::GetVoxelDataSize(std::string key)
{
	NEWMAT::Matrix mat = GetVoxelData(key);
	return mat.Nrows();
}

void FabberRunData::SetVoxelData(std::string key, NEWMAT::Matrix &data)
{
	CheckSize(key, data);
	m_voxel_data[key] = data;
}

void FabberRunData::ClearVoxelData(std::string key)
{
	if (m_voxel_data.count(key) != 0)
	{
		LOG << "FabberRunData::Erasing data " << key << endl;
		m_voxel_data.erase(key);
	}
}

const NEWMAT::Matrix& FabberRunData::GetVoxelData(std::string key)
{
	// Attempt to load data if not already present. Will
	// throw an exception if parameter not specified
	// or file could not be loaded
	//
	// FIXME different exceptions? What about use case where
	// data is optional?
	if (m_voxel_data.count(key) == 0)
	{
#ifdef USE_NEWIMAGE
		string filename = GetStringDefault(key, "");
		LoadVoxelData(filename, key);
#else
		throw DataNotFound(key);
#endif
	}

	return m_voxel_data.find(key)->second;
}

#ifdef USE_NEWIMAGE
#include "newimage/newimageall.h"
using namespace NEWIMAGE;

void DumpVolumeInfo(const volume4D<float>& info, ostream& out = LOG)
{
	Tracer_Plus tr("DumpVolumeInfo");
	LOG << "FabberRunData::Dimensions: x=" << info.xsize() << ", y=" << info.ysize() << ", z=" << info.zsize()
	<< ", vols=" << info.tsize() << endl;
	LOG << "FabberRunData::Voxel size: x=" << info.xdim() << "mm, y=" << info.ydim() << "mm, z=" << info.zdim()
	<< "mm, TR=" << info.tdim() << " sec\n";
	LOG << "FabberRunData::Intents: " << info.intent_code() << ", " << info.intent_param(1) << ", "
	<< info.intent_param(2) << ", " << info.intent_param(3) << endl;
}

void DumpVolumeInfo(const volume<float>& info, ostream& out = LOG)
{
	Tracer_Plus tr("DumpVolumeInfo");
	LOG << "FabberRunData::Dimensions: x=" << info.xsize() << ", y=" << info.ysize() << ", z=" << info.zsize()
	<< ", vols=1" << endl;
	LOG << "FabberRunData::Voxel size: x=" << info.xdim() << "mm, y=" << info.ydim() << "mm, z=" << info.zdim()
	<< "mm, TR=1" << " sec\n";
	LOG << "FabberRunData::Intents: " << info.intent_code() << ", " << info.intent_param(1) << ", "
	<< info.intent_param(2) << ", " << info.intent_param(3) << endl;
}

void FabberRunData::SetVoxelCoordsFromVolume(volume4D<float> first_data)
{
	// First try to get the coords from a mask. This
	// can't be set as voxel data, as it defines what the
	// voxel data is. If present it must be loaded from a file

	// FIXME code duplication
	string mask = GetStringDefault("mask", "");
	if (mask != "")
	{
		LoadVoxelCoordsFromMask(GetString("mask"));
	}
	else
	{
		LOG << "FabberRunData::Setting voxel coordinates from first volume to be loaded" << endl;
		// Create coordinates matrix from the first volume to
		// be loaded
		int nx = first_data.xsize();
		int ny = first_data.ysize();
		int nz = first_data.zsize();
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

		Matrix coords = coordvol.matrix();
		SetVoxelCoords(coords);
	}
	m_have_coords = true;
}

void FabberRunData::LoadVoxelCoordsFromMask(std::string mask_filename)
{
	Tracer_Plus tr("LoadMask");

	LOG_ERR("FabberRunData::Loading mask data from '" + mask_filename << "'" << endl);
	read_volume(m_mask, mask_filename);
	m_mask.binarise(1e-16, m_mask.max() + 1, exclusive);
	DumpVolumeInfo(m_mask);

	// Create coordinates matrix using mask
	LOG << "FabberRunData::Setting voxel coordinates from mask" << endl;

	int nx = m_mask.xsize();
	int ny = m_mask.ysize();
	int nz = m_mask.zsize();
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

	Matrix coords = coordvol.matrix(m_mask);
	SetVoxelCoords(coords);
}

void FabberRunData::LoadVoxelData(std::string filename, std::string key)
{
	Tracer_Plus tr("LoadVoxelData");

	if (filename != "")
	{
		// Load the files.  Note: these functions don't throw errors if file doesn't exist --
		// they just crash.  Hence the detailed logging before we try anything.
		LOG << "FabberRunData::Loading data from '" + filename << "'" << endl;
		volume4D<float> data;
		try
		{
			read_volume4D(data, filename);
		}
		catch (Exception &e)
		{
			throw DataNotFound(key, filename);
		}
		DumpVolumeInfo(data);

		if (!m_have_coords)
		{
			// First data set! Set our co-ords from this
			SetVoxelCoordsFromVolume(data);
		}

		try
		{
			// FIXME m_have_mask would be better
			if (m_mask.xsize() > 0)
			{
				LOG << "     Applying mask to data..." << endl;
				m_voxel_data[key] = data.matrix(m_mask);
			}
			else
			{
				m_voxel_data[key] = data.matrix();
			}
		}
		catch (exception &e)
		{
			LOG << "*** NEWMAT error while thresholding time-series... Most likely a dimension mismatch. ***\n";
			throw e;
		}
	}
	else
	{
		throw DataNotFound(key);
	}
}
#endif

void FabberRunData::GetMainVoxelDataMultiple()
{
	Tracer_Plus tr("GetMainVoxelDataMultiple");

	vector<Matrix> dataSets;
	int n = 1;
	while (true)
	{
		try
		{
			dataSets.push_back(GetVoxelData("data" + stringify(n)));
			n++;
		} catch (DataNotFound &e)
		{
			// No more data sets to combine, carry on with what we've got
			break;
		}
	}

	string order = GetStringDefault("data-order", "interleave");
	int nSets = dataSets.size();
	if (nSets < 1)
	{
		throw Invalid_option(
				"At least one data file is required: --data=<file1> or [--data1=<file1> --data2=<file2> [...]]\n");
	}
	if ((order == "singlefile") && nSets > 1)
	{
		throw Invalid_option("data-order=singlefile but more than one file specified");
	}

	Matrix voxelDataMain;
	if (order == "interleave")
	{
		LOG << "FabberRunData::Combining data into one big matrix by interleaving..." << endl;
		// Interleave - For example if the data sets are A, B, C and each
		// has 3 time points 1, 2, 3 the final time series will be
		// A1B1C1A2B2C2A3B3C3
		int nTimes = dataSets[0].Nrows();
		voxelDataMain.ReSize(nTimes * nSets, dataSets[0].Ncols());
		for (int i = 0; i < nTimes; i++)
		{
			for (int j = 0; j < nSets; j++)
			{
				if (dataSets[j].Nrows() != nTimes)
				{
					// Data sets need same number of time points if they are to be interleaved
					throw Invalid_option("Data sets must all have the same number of time points");
				}
				voxelDataMain.Row(nSets * i + j + 1) = dataSets.at(j).Row(i + 1);
			}
		}
	}
	else if (order == "concatenate")
	{
		LOG << "FabberRunData::Combining data into one big matrix by concatenating..." << endl;
		// Concatentate - For example if the data sets are A, B, C and each
		// has 3 time points 1, 2, 3 the final time series will be
		// A1A2A3B1B2B3C1C2C3
		voxelDataMain = dataSets.at(0);
		for (unsigned j = 1; j < dataSets.size(); j++)
		{
			voxelDataMain &= dataSets.at(j);
		}
	}
	else if (order == "singlefile")
	{
		voxelDataMain = dataSets[0];
	}
	else
	{
		throw Invalid_option("data-order not recognized: " + order);
	}

	m_voxel_data["data"] = voxelDataMain;
	LOG << "FabberRunData::Done loading data, size = " << voxelDataMain.Nrows() << " timepoints by "
			<< voxelDataMain.Ncols() << " voxels" << endl;
}

void FabberRunData::SaveVoxelData(std::string filename, NEWMAT::Matrix &data, int nifti_intent_code)
{
	int data_size = data.Nrows();

	if (!m_save_files)
	{
		LOG << "FabberRunData::Saving to memory: " << filename << endl;
		// FIXME what should we do with NIFTI_INTENT_CODE?
		SetVoxelData(filename, data);
	}
	else
	{
#ifdef USE_NEWIMAGE
		LOG << "FabberRunData::Saving to nifti: " << filename << endl;
		volume4D<float> output(m_size[0], m_size[1], m_size[2], data_size);
		if (m_mask.xsize() > 0)
		{
			output.setmatrix(data, m_mask);
		}
		else
		{
			output.setmatrix(data);
		}
		output.set_intent(nifti_intent_code, 0, 0, 0);
		output.setDisplayMaximumMinimum(output.max(), output.min());
		string output_dir = GetStringDefault("output", ".");
		save_volume4D(output, output_dir + "/" + filename);
#else
		throw Invalid_option("Asked to save data to file, but file I/O via NEWIMAGE not supported in this version")
#endif
	}
}

