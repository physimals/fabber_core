/*  dataset.h - Data-loading class for FABBER

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "fabber_io.h"
#include "easylog.h"

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
 *
 * Mostly self-explanatory, apart from:
 *
 * OPT_IMAGE is 3D single-valued voxel data, e.g. image prior
 * OPT_TIMESERIES is 4D voxel data, e.g. main/supplemental data.
 * OPT_MVN is MVN voxel data
 * OPT_MATRIX is a matrix file (either VEST or ASCII), usually used for small matrices
 * OPT_FILE is some other type of file
 */
enum OptionType
{
	OPT_BOOL, OPT_STR, OPT_INT, OPT_FLOAT, OPT_FILE, OPT_IMAGE, OPT_TIMESERIES, OPT_MVN, OPT_MATRIX
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
	virtual void Progress(int voxel, int nVoxels)
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
	virtual void Progress(int voxel, int nVoxels);

private:
	int m_last;
};

/**
 * ProgressCheck which calls a C function pointer
 *
 * Used for the C API
 */
class CallbackProgressCheck: public ProgressCheck
{
public:
	CallbackProgressCheck(void (*cb)(int, int)) :
			m_cb(cb)
	{
	}
	virtual void Progress(int voxel, int nVoxels) {m_cb(voxel, nVoxels);}

private:
	void (*m_cb)(int, int);
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
class FabberRunData: public Loggable
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

	/**
	 * Constructor
	 *
	 * @param io Instance responsible for loading/saving voxel data.
	 *           This will not be copied or freed. The caller is responsible
	 *           for freeing it after use.
	 */
	FabberRunData(FabberIo *io, bool compat_options = true);
	~FabberRunData();

	/**
	 * Run fabber
	 */
	void Run(ProgressCheck *check = 0);

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
	 * 'Unset' an option
	 *
	 * If the option is normally treated as a bool its
	 * value will now be 'false'
	 */
	void Unset(const std::string key);

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
	 * Get the output directory for this run.
	 *
	 * This is derived from the value of the 'output' parameter.
	 * However if this directory already exists and the 'overwrite'
	 * boolean parameter is not set, then we append '+' to the
	 * directory name until we obtain a unique directory, which
	 * will then be created. If more than 50 such directories are
	 * tried, give up and throw an exception
	 *
	 * If 'overwrite' is specified, we just use the directory specified
	 * in the 'output' parameter regardless of whether it exists or
	 * not. If it cannot be written to an exception is thrown.
	 *
	 * If 'output' is not set, "." is returned so output is written
	 * to the working directory.
	 */
	std::string GetOutputDir();

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
	const NEWMAT::Matrix& GetVoxelCoords() const;

	/**
	 * Get named voxel data
	 *
	 * GetVoxelData will check recursively for a string option with this key
	 * until it can resolve the key no further. The last non-empty key will
	 * then be used to request data from the FabberIo instance.
	 *
	 * For example, if GetVoxelData("data") is called, we first check for a
	 * string option with key "data". If this option is set to "fmri_data",
	 * then we check for a string option with key "fmri_data", and so on.
	 * If no option with key "fmri_data" is found, we ask our FabberIo instance
	 * to get the data with key "fmri_data", which it might, for example,
	 * load from a file fmri_data.nii.gz.
	 *
	 * This sounds complicated, but actually simplifies situations such as
	 * restarting runs from a file, where the restart might come from memory
	 * saved data, or an external file.
	 *
	 * @param key Name identifying the voxel data required.
	 * @return an NxM matrix where each column contains the data for a single
	 *         voxel. The rows may contain the time series of data for that voxel
	 *         however they might be used for other purposes, e.g. storing the mean
	 *         of each parameter for that voxel.
	 * @throw DataNotFound If no voxel data matching key is found and no data
	 *                     could be loaded
	 */
	const NEWMAT::Matrix& GetVoxelData(std::string key, bool allowFile = false);

	/**
	 * Get the number of data values associated with each voxel for the named data
	 *
	 * @param key Name identifying the voxel data
	 * @return The number of data values, e.g. 1 for a simple image, or for a time
	 *         series the number of time slices, etc.
	 */
	int GetVoxelDataSize(std::string key);

	/**
	 * Get the main voxel data
	 *
	 * This will initially call GetVoxelData with the key "data", however if
	 * this is not found it may attempt to load data from multiple files
	 * "data1", "data2", etc, using the data-order option.
	 *
	 * @return an NxT matrix where each column contains the data for a single
	 *         voxel. The rows contain the time series of data for that voxel
	 */
	const NEWMAT::Matrix& GetMainVoxelData();

	/**
	 * Get the voxel supplementary data
	 *
	 * This is equivalent to calling GetVoxelData("suppdata"), except that if
	 * no supplemental data is found this method returns an empty matrix rather
	 * than throwing an exception.
	 *
	 * @return an NxT matrix where each column contains the supplementary data for a single
	 *         voxel. The rows contain the time series of data for that voxel
	 */
	const NEWMAT::Matrix& GetVoxelSuppData();

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
			m_progress->Progress(voxel, nVoxels);
	}

	/**
	 * Get the I/O instance being used
	 */
	FabberIo* GetIo()
	{
		return m_io;
	}

	/**
	 * Send list of all parameters to the logfile
	 */
	void LogParams();

	/**
	 * Friend function to allow summary of Parameters to be streamed
	 * using the << operator.
	 */
	friend ostream& operator<<(ostream& out, const FabberRunData& opts);

	// Following methods present for compatibility only
#ifdef DEPRECATED
	/** @deprecated Use GetString instead */
	std::string Read(std::string key);

	/** @deprecated Use GetString instead */
	std::string Read(std::string key, std::string msg);

	/** @deprecated Use GetStringDefault instead */
	std::string ReadWithDefault(std::string key, std::string def);

	/** @deprecated Use GetBool instead */
	bool ReadBool(std::string key);

	void ParseOldStyleParamFile(const string filename);
#endif

private:
	void init(bool compat_options);
	void AddKeyEqualsValue(const std::string key, bool trim_comments = false);
	const NEWMAT::Matrix &GetMainVoxelDataMultiple();

	/** IO instance - not owned and will not be freed */
	FabberIo *m_io;

	/** Optional progress checker, could be NULL - not owned and will not be freed */
	ProgressCheck *m_progress;

	/**
	 * Empty matrix
	 *
	 * Used for returning SUPPDATA which is optional and methods expect to receive
	 * and empty matrix when it is not supplied
	 */
	NEWMAT::Matrix m_empty;

	/**
	 * Stores main voxel data when supplied in concatenated/interleaved form
	 */
	NEWMAT::Matrix m_mainDataMultiple;

	/**
	 * Options as key/value pairs
	 */
	std::map<std::string, std::string> m_params;

	std::string m_outdir;

	EasyLog m_default_log;
};

/**
 * Exception base class
 */
class FabberError: public std::runtime_error
{
public:
	FabberError(std::string msg) :
			std::runtime_error(msg), m_msg(msg)
	{
	}

	virtual ~FabberError() throw() {};
	virtual const char* what() const throw() {return m_msg.c_str();}

	std::string m_msg;
};

/**
 * Thrown for errors that are not the direct
 * result of bad user input, i.e. bugs in fabber
 */
class FabberInternalError: public FabberError
{
public:
	FabberInternalError(std::string msg) :
		FabberError(msg)
	{
		m_msg = "Internal error in Fabber: " + msg;
	}
};

/**
 * Thrown for errors that are directly caused
 * by invalid user input
 */
class FabberRunDataError: public FabberError
{
public:
	FabberRunDataError(std::string msg) :
		FabberError(msg)
	{
	}
};

/**
 * Thrown when the value given to an option is not valid,
 * e.g. --num_mc_steps=fred (expected a number)
 */
class InvalidOptionValue: public FabberRunDataError
{
public:
	InvalidOptionValue(std::string key, std::string value, std::string reason = "") :
		FabberRunDataError(key)
	{
		m_msg = "Invalid value given for option: " + key + "=" + value + " (" + reason + ")";
	}
};

/**
 * Thrown when a mandatory option or parameter is not given
 */
class MandatoryOptionMissing: public FabberRunDataError
{
public:
	MandatoryOptionMissing(std::string key) :
		FabberRunDataError(key)
	{
	  m_msg = "No value given for mandatory option: " + key;
	}
};

/**
 * Thrown when voxel data is requested but cannot be found
 */
class DataNotFound: public FabberRunDataError
{
public:
	DataNotFound(std::string key, std::string reason="") :
		FabberRunDataError(key)
	{
		m_msg = "Voxel data not found: " + key + " (" + reason + ")";
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
		throw InvalidOptionValue("", s, "Failed to convert to required type");
	return x;
}

#ifdef DEPRECATED
typedef class FabberRunData ArgsType;
typedef class FabberRunData EasyOptions;
typedef class FabberError Invalid_option;
#endif
