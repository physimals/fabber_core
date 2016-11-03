/*  dataset.cc - Data-loading class for FABBER

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "dataset.h"

#include "fabber_version.h"
#include "easylog.h"
#include "setup.h"
#include "fwdmodel.h"
#include "inference.h"
#include "fabber_io.h"

#include "utils/tracer_plus.h"
#include "newmat.h"

#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using Utilities::Tracer_Plus;

std::ostream& operator<<(std::ostream& out, const OptionType value)
{
	const char* s = 0;
	switch (value)
	{
	case OPT_BOOL:
		s = "BOOL";
		break;
	case OPT_STR:
		s = "STRING";
		break;
	case OPT_INT:
		s = "INTEGER";
		break;
	case OPT_FLOAT:
		s = "FLOAT";
		break;
	case OPT_FILE:
		s = "FILENAME";
		break;
	case OPT_MATRIX:
		s = "MATRIXFILE";
		break;
	default:
		s = "UNKNOWN";
	}

	return out << s;
}

std::ostream& operator<<(std::ostream& out, const OptionSpec &value)
{
	return out << "--" << value.name << " [" << value.type << "," << (value.optional ? "NOT REQUIRED" : "REQUIRED")
			<< "," << ((value.def == "") ? "NO DEFAULT" : "DEFAULT=" + value.def) << "]" << endl << "        "
			<< value.description << endl;
}

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

string FabberRunData::GetVersion()
{
	stringstream v;
	v << V_MAJ << "." << V_MIN << "." << V_PAT << V_FL;
	return v.str();
}

string FabberRunData::GetRevision()
{
#ifdef GIT_SHA1
	return GIT_SHA1;
#else
	return "unknown";
#endif
}

string FabberRunData::GetDate()
{
#ifdef GIT_DATE
	return GIT_DATE;
#else
	return "unknown";
#endif
}

FabberRunData::FabberRunData(FabberIo *io)
 : m_save_files(false), m_have_coords(false), m_nvoxels(-1), m_progress(0), m_io(io)
{
	FabberSetup::SetupDefaults();
}

FabberRunData::FabberRunData() :
		m_save_files(false), m_have_coords(false), m_nvoxels(-1), m_progress(0), m_io(0)
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

	LOG << "FabberRunData::FABBER release v" << GetVersion() << endl;
	LOG << "FabberRunData::Revision " << GetRevision() << endl;
	LOG << "FabberRunData::Last commit: " << GetDate() << endl;

	time_t startTime;
	time(&startTime);
	LOG << "FabberRunData::Start time: " << ctime(&startTime);

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

	// FIXME this is a hack but seems to be expected that the command line
	// tool will output parameter names to a file. Really should be an option!
	if (m_save_files)
	{
		ofstream paramFile((EasyLog::GetOutputDirectory() + "/paramnames.txt").c_str());
		vector<string> paramNames;
		fwd_model->NameParams(paramNames);
		for (unsigned i = 0; i < paramNames.size(); i++)
		{
			paramFile << paramNames[i] << endl;
		}
		paramFile.close();
	}

	LOG << "FabberRunData::All done." << endl;

	time_t endTime;
	time(&endTime);
	LOG << "FabberRunData::Start time: " << ctime(&startTime); // Bizarrely, ctime() ends with a \n.
	LOG << "FabberRunData::End time: " << ctime(&endTime);
	LOG << "FabberRunData::Duration: " << endTime - startTime << " seconds." << endl;
}

static string trim(string const& str)
{
	if (str.empty())
		return str;

	size_t firstScan = str.find_first_not_of(' ');
	size_t first = firstScan == string::npos ? str.length() : firstScan;
	size_t last = str.find_last_not_of(' ');
	return str.substr(first, last - first + 1);
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
	else if (eqPos != (exp.npos))
	{
		string::size_type end = exp.npos;
		if (trim_comments)
			end = exp.find("#");
		string value = trim(exp.substr(eqPos + 1, end - (eqPos + 1)));
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

int FabberRunData::GetInt(const string key)
{
	string val = GetString(key, "Mandatory option not specified: " + key);
	try
	{
		return convertTo<int>(val);
	} catch (invalid_argument&)
	{
		throw Invalid_option(key + " must be an integer");
	}
}

double FabberRunData::GetDouble(const string key)
{
	string val = GetString(key, "Mandatory option not specified: " + key);
	try
	{
		return convertTo<double>(val);
	} catch (invalid_argument&)
	{
		throw Invalid_option(key + " must be a number (was: " + val + ")");
	}
}

int FabberRunData::GetIntDefault(const string key, int def)
{
	if (m_params.count(key) == 0)
		return def;
	else
		return GetInt(key);
}

double FabberRunData::GetDoubleDefault(const string key, double def)
{
	if (m_params.count(key) == 0)
		return def;
	else
		return GetDouble(key);
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

	// FIXME we assume coords will not be negative
	m_size.resize(3);
	m_size[0] = coords.Row(1).Maximum() - coords.Row(1).Minimum() + 1;
	m_size[1] = coords.Row(2).Maximum() - coords.Row(2).Minimum() + 1;
	m_size[2] = coords.Row(3).Maximum() - coords.Row(3).Minimum() + 1;

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
		string filename = GetStringDefault(key, "");
		if (m_io && (filename != ""))
		{
			// If this is the first data to be loaded, load any
			// mask first
			if (!m_have_coords)
			{
				string mask_filename = GetStringDefault("mask", "");
				if (mask_filename != "") m_io->LoadMask(mask_filename);
			}
			m_voxel_data[key] = m_io->LoadVoxelData(filename);

			// Set the coords last as they may have come either
			// from the mask or from the data
			if (!m_have_coords) {
				Matrix coords = m_io->GetVoxelCoords();
				SetVoxelCoords(coords);
			}
		}
		else
		{
			throw DataNotFound(key);
		}
	}

	return m_voxel_data.find(key)->second;
}

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

void FabberRunData::SaveVoxelData(std::string filename, NEWMAT::Matrix &data, VoxelDataType data_type)
{
	if (m_io && m_save_files)
	{
		m_io->SaveVoxelData(data, m_size, filename, data_type);
	}
	else
	{
		LOG << "FabberRunData::Saving to memory: " << filename << endl;
		// FIXME what should we do with data_type?
		SetVoxelData(filename, data);
	}
}

