/*  dataset.h - Data-loading class for FABBER

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "fabber_io.h"

#include "newmat.h"

#include <boost/shared_ptr.hpp>

#include <stdexcept>
#include <vector>
#include <map>
#include <sstream>

/** Include deprecated compatibility methods */
#define DEPRECATED 7

/**
 * Option types
 */
enum OptionType
{
	OPT_BOOL, OPT_STR, OPT_INT, OPT_FLOAT, OPT_FILE, OPT_MATRIX
};

/**
 * Option required
 */
enum OptionReq
{
	OPT_REQ = 0, OPT_NONREQ = 1,
};

std::ostream& operator<<(std::ostream& out, const OptionType value);

/**
 * Describes a runtime option for a model, inference method, etc
 */
struct OptionSpec
{
	/** Option name as used on command line --name */
	std::string name;
	/** Option type - string or boolean */
	OptionType type;
	/** Human readable description */
	std::string description;
	/** Whether option needs to be specified */
	OptionReq optional;
	/** Default value if there is one */
	std::string def;
};

std::ostream& operator<<(std::ostream& out, const OptionSpec &value);

/**
 * Functor which is be called to monitor progress.
 *
 * The default does nothing.
 */
class ProgressCheck
{
public:
	virtual void operator()(int voxel, int nVoxels)
	{
	}
};

/**
 * ProgressCheck which prints percentages in multiples of 10 to stdout
 *
 * This is basically an example of how to use it
 */
class PercentProgressCheck: public ProgressCheck
{
public:
	PercentProgressCheck() :
			m_last(-1)
	{
	}
	virtual void operator()(int voxel, int nVoxels);

private:
	int m_last;
};

/**
 * Encapsulates all the input and output data associated with a fabber run
 *
 * Run data is of the following types:
 *   - String parameters with a key and value
 *   - Boolean parameters with a key which are either present (true) or not (false)
 *   - Voxel data, which is a set of N values associated with each voxel
 *
 * In a command line program, run data is loaded during parsing of the command
 * line options in the Parse() or ParseParamFile methods.
 *
 * In a library program, it should be provided explicitly by calling the various
 * Set() methods.
 *
 * The class will try to keep all it's voxel data consistent. Once any data is
 * obtained, subsequent per-voxel data sets must be of the correct size.
 */
class FabberRunData
{
public:

	/**
	 * Get current version string
	 *
	 * @return A version string in form major.minor.patch[_flag], e.g. 1.2.3_rc1
	 */
	static std::string GetVersion();

	/**
	 * Get SHA1 hash of current Git revision, if available
	 *
	 * @return an SHA1 hash, or 'unknown' if not available at build time
	 */
	static std::string GetRevision();

	/**
	 * Get date of last commit, if available
	 *
	 * @return date string, or 'unknown' if not available at build time
	 */
	static std::string GetDate();

	/**
	 * Get general Fabber option descriptions
	 */
	static void GetOptions(std::vector<OptionSpec> &opts);

	FabberRunData(FabberIo *io);
	FabberRunData();
	~FabberRunData();

	/**
	 * Run fabber
	 */
	void Run();

	/**
	 * Parse command line arguments into run data
	 *
	 * Accepted forms:
	 * --key=value -> args[key] == value
	 * --key value -> args[key] == value
	 * --key       -> args[key] == ""
	 * For the last one, use args.count(key) as a boolean flag.  Do not look at
	 * args[key] because that will create the empty key for you!
	 * Also, the command name is stored in args[""].

	 * Extensions:
	 * - use a multimap so duplicate options are okay
	 * - accept short-form options (e.g. -k)
	 */
	void Parse(int argc, char** argv);

	/**
	 * Parse options from a .fab file
	 *
	 * This file has a simple format:
	 *
	 * # This is a comment
	 * option=value
	 * bool-option
	 */
	void ParseParamFile(const std::string file);

	/**
	 * Set whether to save files or not
	 *
	 * @param save if true, voxel data which is saved by
	 * inference methods will be written to corresponding
	 * files. Otherwise it will be kept in the run data
	 * where it can be accessed using GetVoxelData()
	 */
	void SetSaveFiles(bool save)
	{
		m_save_files = save;
	}

	/**
	 * Set string option.
	 *
	 * Will overwrite any boolean option of the same name
	 *
	 * @param key Name of the string option
	 * @param value It's value
	 */
	void Set(const std::string key, const std::string value);

	/**
	 * Set string option from a number.
	 *
	 * Will overwrite any boolean option of the same name
	 *
	 * @param key Name of the string option
	 * @param value which will be converted to string
	 */
	void Set(const std::string key, double value);

	/**
	 * Set boolean option.
	 *
	 * Will overwrite any string option of the same name
	 *
	 * Note we do not overload the Set() method because
	 * C++ will then misinterpret any call to Set which
	 * uses a C-string as the value.
	 *
	 * @param key Name of the boolean option
	 * @param value true or false
	 */
	void SetBool(const std::string key, bool value = true);

	/**
	 * Get string option.
	 *
	 * @param key Name of the string option
	 * @throw if option is missing or is a boolean
	 */
	std::string GetString(const std::string key);

	/**
	 * Get string option.
	 *
	 * @param key Name of the string option
	 * @param msg Custom error message if not found
	 * @throw if option is missing or is a boolean
	 */
	std::string GetString(const std::string key, const std::string msg);

	/**
	 * Get string option with default if not found.
	 *
	 * @param key Name of the string option
	 * @param def Default value to return if not found
	 */
	std::string GetStringDefault(const std::string key, const std::string def);

	/**
	 * Get boolean option
	 *
	 * @param key Name of the boolean option
	 * @throw if string option specified
	 */
	bool GetBool(const std::string key);

	/**
	 * Get an integer option
	 *
	 * @param key Name of the option
	 * @throw if option not specified, or not an integer
	 */
	int GetInt(const std::string key);

	/**
	 * Get an integer option, returning default if not specified
	 *
	 * @param key Name of the option
	 * @throw if option specified, but not a valid integer
	 */
	int GetIntDefault(const std::string key, int def);

	/**
	 * Get an double option
	 *
	 * @param key Name of the option
	 * @throw if option not specified, or not a double
	 */
	double GetDouble(const std::string key);

	/**
	 * Get an double option, returning default if not specified
	 *
	 * @param key Name of the option
	 * @throw if option specified, but not a valid number
	 */
	double GetDoubleDefault(const std::string key, double def);

	/**
	 * Check that all key/value parameters have been consumed
	 *
	 * Used by the command line tool to detect if options
	 * have been specified which have never been used by
	 * the program. This usually suggests a user error
	 */
	void CheckEmpty();

	//void LoadVest(std::string filename, std::string key);

	/**
	 * Save the specified voxel data
	 *
	 * If SetSaveFiles has been set to true, the specified data will be written
	 * out to a file but not kept in memory. This is the default when
	 * command line options have been used to initialize the rundata using Parse
	 *
	 * If SetSaveFiles is fales, the data will be saved internally
	 * and can be retrieved using GetVoxelData. The filename specified is the
	 * data key value.
	 *
	 * @param filename Filename to write to, or the key name for the data
	 * @param data Data as a matrix in which each column is a voxel, and
	 *        rows contain a series of data values for that voxel
	 */
	void SaveVoxelData(std::string filename, NEWMAT::Matrix &coords, VoxelDataType data_type = VDT_SCALAR);

	/**
	 * Get the voxel co-ordinates
	 *
	 * @return an Nx3 matrix where each column is a voxel and the rows
	 *         are the xyz co-ords of the voxel. The co-ordinates are
	 *         grid positions (integers), not physical co-ordiantes (mm)
	 */
	const NEWMAT::Matrix& GetVoxelCoords() const
	{
		return m_voxelCoords;
	}

	/**
	 * Set the voxel co-ordinates
	 *
	 * @param coords an Nx3 matrix where each column is a voxel and the rows
	 *         are the xyz co-ords of the voxel. The co-ordinates are
	 *         grid positions (integers), not physical co-ordiantes (mm)
	 */
	void SetVoxelCoords(NEWMAT::Matrix &coords);

	/**
	 * Get named voxel data
	 *
	 * @param key Name identifying the voxel data required. If no voxel data
	 *            matching key is found, will attempt to load data from
	 *            a NIFTI file using the string parameter matching key.
	 * @return an NxM matrix where each column contains the data for a single
	 *         voxel. The rows may contain the time series of data for that voxel
	 *         however they might be used for other purposes, e.g. storing the mean
	 *         of each parameter for that voxel.
	 * @throw DataNotFound If no voxel data matching key is found and no data
	 *                     could be loaded
	 */
	const NEWMAT::Matrix& GetVoxelData(std::string key);

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
	void SetVoxelData(std::string key, NEWMAT::Matrix &data);

	/**
	 * Remove previously set voxel data
	 *
	 * If no such data exists, does nothing
	 *
	 * @key Name identifying the voxel data to clear
	 */
	void ClearVoxelData(std::string key);

	/**
	 * Get the number of data values associated with each voxel for the named data
	 *
	 * @param key Name identifying the voxel data
	 * @return The number of data values, e.g. 1 for a simple image, or for a time
	 *         series the number of time slices, etc.
	 */
	int GetVoxelDataSize(std::string key);

	/**
	 * Get the number of voxels in the data
	 */
	int GetNVoxels() const
	{
		return m_nvoxels;
	}

	/**
	 * Get the size of the volume grid, i.e. the number of x, y and z values
	 * along each axis
	 */
	std::vector<int> GetVolumeSize() const
	{
		return m_size;
	}

	/**
	 * Get the size of each voxel in mm in the x, y and z directions;
	 */
	std::vector<float> GetVoxelDims() const
	{
		return m_dims;
	}

	/**
	 * Get the main voxel data
	 *
	 * @return an NxT matrix where each column contains the data for a single
	 *         voxel. The rows contain the time series of data for that voxel
	 */
	const NEWMAT::Matrix& GetMainVoxelData();

	/**
	 * Set the main voxel data
	 *
	 * @param coords an NxT matrix where each column contains the data for a single
	 *         voxel. The rows contain the time series of data for that voxel
	 */
	void SetMainVoxelData(NEWMAT::Matrix &data);

	/**
	 * Get the voxel supplementary data
	 *
	 * @return an NxT matrix where each column contains the supplementary data for a single
	 *         voxel. The rows contain the time series of data for that voxel
	 */
	const NEWMAT::Matrix& GetVoxelSuppData();

	/**
	 * Set the voxel supplementary data
	 *
	 * @param coords an NxT matrix where each column contains supplementary data for a single
	 *         voxel. The rows contain the time series of data for that voxel
	 */
	void SetVoxelSuppData(NEWMAT::Matrix &data);

	/**
	 * Pass a functor which will be called periodically to allow progress to
	 * be reported
	 *
	 * The only guarantees are that the functor will be called once at the
	 * start and once at the end.
	 */
	void SetProgressCheck(ProgressCheck *check)
	{
		m_progress = check;
	}

	/**
	 * Report progress
	 *
	 * InferenceMethods call this to report how many voxels
	 * have been completed. It will call the user-supplied
	 * progress checking functor
	 *
	 * @param voxel Current voxel
	 * @param nVoxels Total number of voxels
	 */
	void Progress(int voxel, int nVoxels)
	{
		if (m_progress)
			(*m_progress)(voxel, nVoxels);
	}

	/**
	 * Friend function to allow summary of Parameters to be streamed
	 * using the << operator.
	 */
	friend ostream& operator<<(ostream& out, const FabberRunData& opts);

	// Following methods present for compatibility only
#ifdef DEPRECATED
	/** @deprecated Use GetString instead */
	std::string Read(std::string key)
	{
		return GetString(key);
	}
	/** @deprecated Use GetString instead */
	std::string Read(std::string key, std::string msg)
	{
		return GetString(key);
	}
	/** @deprecated Use GetStringDefault instead */
	std::string ReadWithDefault(std::string key, std::string def)
	{
		return GetStringDefault(key, def);
	}
	/** @deprecated Use GetBool instead */
	bool ReadBool(std::string key)
	{
		return GetBool(key);
	}
	void ParseOldStyleParamFile(const string filename);
#endif

	/**
	 * Send list of all parameters to the logfile
	 */
	void LogParams();
private:
	void init();
	void AddKeyEqualsValue(const std::string key, bool trim_comments = false);
	void CheckSize(std::string key, NEWMAT::Matrix &mat);
	void GetMainVoxelDataMultiple();

	FabberIo *m_io;
	NEWMAT::Matrix m_voxelCoords;
	bool m_save_files;
	bool m_have_coords;
	int m_nvoxels;
	NEWMAT::Matrix m_empty;
	std::map<std::string, std::string> m_params;
	std::vector<int> m_size;
	std::vector<float> m_dims;
	std::map<std::string, NEWMAT::Matrix> m_voxel_data;
	std::map<std::string, NEWMAT::Matrix> m_misc_data;

	ProgressCheck *m_progress;
};

#ifdef DEPRECATED
typedef class FabberRunData ArgsType;
typedef class FabberRunData EasyOptions;
#endif

/**
 * Thrown when the value given to an option is not valid,
 * e.g. --num_mc_steps=fred (expected a number)
 */
class Invalid_option: public std::runtime_error
{
public:
	Invalid_option(std::string msg) :
			std::runtime_error(msg)
	{
	}
};

/**
 * Thrown when a mandatory option or parameter is not given
 */
class MandatoryParameterMissing: public Invalid_option
{
public:
	MandatoryParameterMissing(std::string msg) :
			Invalid_option(msg)
	{
	}
};

/**
 * Thrown when data is requested but cannot be found
 */
class DataNotFound: public Invalid_option
{
public:
	DataNotFound(std::string key, std::string filename) :
			Invalid_option(key + " (filename: " + filename + ")")
	{
	}
	DataNotFound(std::string key) :
			Invalid_option(key)
	{
	}
};

// Convert a string into almost anything, using operator>>. See the C++ FAQ-Lite:
// http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.3
template<typename T> inline T convertTo(const std::string& s)
{
	T x;
	istringstream i(s);
	char c;
	if (!(i >> x) || (i.get(c)))
		throw Invalid_option("convertTo failed... couldn't convert string '" + s + "' into a value.\n");
	return x;
}
