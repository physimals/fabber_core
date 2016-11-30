/*  easylog.h - a fairly minimal logging-to-file implementation

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include "assert.h"

#define PRINTNOTE fprintf(stderr, "Note: %s line %d\n", __FILE__, __LINE__);

/**
 * Use LOG just like you'd use cout.  (i.e. LOG << your_data << endl;)
 *
 * StartLog() is used to tell the logger where to put the output - any log
 * calls made before this are stored in a temporary stringstream and then
 * flushed to the log once StartLog is called.
 *
 * If you want to dump variables directly to this stream, just implement
 * std::ostream& operator<<(std::ostream& stream, const your_type& data);
 */
#define LOG (EasyLog::CurrentLog().LogStream())

#define LOG_ERR(x) ((void)((LOG<<x)))

#define WARN_ONCE(x) ((void)(EasyLog::CurrentLog().WarnOnce(x)))

#define WARN_ALWAYS(x) ((void)(EasyLog::CurrentLog().WarnAlways(x)))

/**
 * Sends logging information to an output stream
 *
 * The stream might be a file, a string buffer, or
 * stdout/stderr
 */
class EasyLog
{
public:

	/**
	 * Return pointer to current log.
	 */
	static EasyLog& CurrentLog();

	/**
	 * Log output to an existing stream
	 *
	 * This might be std::cout/cerr or something else, e.g. a
	 * stringstream.
	 */
	EasyLog();

	/**
	 * Log output to a file
	 *
	 * The file will be called 'logfile'
	 *
	 * @param basename Name of a directory in which to put the logfile.
	 * @param overwrite If true, will replace any existing log file in the
	 *                  specified directory. If false, and the specified
	 *                  directory already exists, a new directory will
	 *                  be created with a '+' added to basename, e.g.
	 *                  'out', 'out+', 'out++', etc...
	 */
	void StartLog(const std::string& basename, bool overwrite, bool link_to_latest=false);

	/**
	 * Log output to an existing stream
	 *
	 * This might be std::cout/cerr or something else, e.g. a
	 * stringstream.
	 */
	void StartLog(std::ostream& s);

	/**
	 * Get the output directory for the log file, or
	 * an empty string if not logging to a file.
	 */
	const std::string& GetOutputDirectory();

	/**
	 * Stop logging.
	 *
	 * @param gzip If true and we are logging to a file, gzip the logfile
	 */
	void StopLog(bool gzip = false);

	/**
	 * @return true if StartLog has been called
	 */
	bool LogStarted();

	/**
	 * Get the logging stream. Easier to use the LOG macro defined above
	 */
	std::ostream& LogStream();

	/**
	 * Issue a warning
	 *
	 * The warning will appear once in the log.
	 * If the same warning occurs again, it will not
	 * be repeated, apart from if ReissueWarnings is
	 * called
	 */
	void WarnOnce(const std::string& text);

	/**
	 * Issue a warning
	 *
	 * The warning will appear in the log each time it is issued.
	 */
	void WarnAlways(const std::string& text);

	/**
	 * Resend all warnings recorded so far to the log
	 */
	void ReissueWarnings();

private:
	std::ostream* filestream;
	std::stringstream templog;
	std::string outDir;
	std::map<std::string, int> warnCount;

	// Thread-local loggers
#ifdef _WIN32
	static map<DWORD, EasyLog> s_logs;
#else
  #ifdef USE_PTHREADS
	static map<pthread_t , EasyLog> s_logs;
  #else
	static EasyLog s_log;
  #endif
#endif
};

// Other useful functions:

// Convert almost anything to a string
#include <sstream>
#include <stdexcept>
template<typename type>
inline std::string stringify(type from)
{
	std::ostringstream s;
	if (!(s << from))
		throw std::logic_error("Stringify failed");
	return s.str();
}

// output a vector<int>
inline std::ostream& operator<<(std::ostream& out, std::vector<int> x)
{
	out << "[ ";
	for (unsigned i = 0; i < x.size(); i++)
		out << x[i] << " ";
	out << "]";
	return out;
}
