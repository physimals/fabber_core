/*  easylog.h - a fairly minimal logging-to-file implementation

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include "assert.h"

#include "utils/tracer_plus.h"

using namespace std;
using namespace Utilities;

#define PRINTNOTE fprintf(stderr, "Note: %s line %d\n", __FILE__, __LINE__);

/**
 * Use LOG just like you'd use cout.  (i.e. LOG << your_data << endl;)
 *
 * StartLog() is used to tell the logger where to put the output - any log
 * calls made before this are stored in a temporary stringstream and then
 * flushed to the log once StartLog is called.
 *
 * If you want to dump variables directly to this stream, just implement
 * ostream& operator<<(ostream& stream, const your_type& data);
 */
#define LOG (*EasyLog::CurrentLog())

#define LOG_ERR(x) ((void)((LOG<<x)))

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
	 * Returns a pointer to the current logging stream
	 *
	 * The result can be used with the << operator, e.g.
	 *   *CurrentLog() << "Hello World" << endl;
	 */
	static ostream* CurrentLog()
	{
		if (filestream == NULL)
		{
			return &templog;
		}
		else
		{
			return filestream;
		}
	}

	/**
	 * Get the output directory for the log file, or
	 * an empty string if not logging to a file.
	 */
	static const string& GetOutputDirectory()
	{
		assert(filestream != NULL);
		return outDir;
	}

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
	static void StartLog(const string& basename, bool overwrite);

	/**
	 * Log output to an existing stream
	 *
	 * This might be std::cout/cerr or something else, e.g. a
	 * stringstream.
	 */
	static void StartLog(ostream& s);

	/**
	 * Stop logging.
	 *
	 * @param gzip If true and we are logging to a file, gzip the logfile
	 */
	static void StopLog(bool gzip = false);

	/**
	 * @return true if StartLog has been called
	 */
	static bool LogStarted()
	{
		return filestream != NULL;
	}

private:
	static stringstream templog;
	static ostream* filestream;
	static string outDir;
};

// Other useful functions:

// Convert almost anything to a string
#include <sstream>
#include <stdexcept>
template<typename type>
inline string stringify(type from)
{
	ostringstream s;
	if (!(s << from))
		throw logic_error("Stringify failed");
	return s.str();
}

// output a vector<int>
#include <vector>
inline ostream& operator<<(ostream& out, vector<int> x)
{
	out << "[ ";
	for (unsigned i = 0; i < x.size(); i++)
		out << x[i] << " ";
	out << "]";
	return out;
}

#include <map>

/**
 * Allows code to issue warnings which can be recorded
 * and repeated at the end so they do not get missed
 */
class Warning
{
public:
	/**
	 * Issue a warning
	 *
	 * The warning will appear once in the current log.
	 * If the same warning occurs again, it will not
	 * be repeated.
	 */
	static void IssueOnce(const string& text);

	/**
	 * Issue a warning
	 *
	 * The warning will appear in the current log
	 * each time it is issued.
	 */
	static void IssueAlways(const string& text);

	/**
	 * Resend all warnings recorded so far to the current log
	 */
	static void ReissueAll();
private:
	static map<string, int> issueCount;
};
