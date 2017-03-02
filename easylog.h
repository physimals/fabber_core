/*  easylog.h - a fairly minimal logging-to-file implementation

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include <assert.h>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * Use LOG just like you'd use cout.  (i.e. LOG << your_data << endl;)
 * Your class must be derived from Loggable.
 *
 * If the log member variable m_log is not set, log output goes to stderr
 */
#define LOG ((m_log == 0) ? std::cerr : m_log->LogStream())
#define LOG_ERR(x) ((void)((LOG << x)))

#define WARN_ONCE(x) \
    if (m_log)       \
    m_log->WarnOnce(x)
#define WARN_ALWAYS(x) \
    if (m_log)         \
    m_log->WarnAlways(x)

/**
 * Sends logging information to an output stream
 *
 * StartLog() is used to tell the logger where to put the output - any log
 * calls made before this are stored in a temporary stringstream and then
 * flushed to the log once StartLog is called.
 *
 * If you want to dump variables directly to this stream, just implement
 * std::ostream& operator<<(std::ostream& stream, const your_type& data);
 */
class EasyLog {
public:
    /**
	 * Log output to an existing stream
	 *
	 * This might be std::cout/cerr or something else, e.g. a
	 * stringstream.
	 */
    EasyLog();
    ~EasyLog();

    /**
	 * Log output to a file
	 *
	 * The file will be called 'logfile'
	 *
	 * @param outDir Name of a directory in which to put the logfile.
	 */
    void StartLog(const std::string& outDir);

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
    std::ostream* stream;
    std::stringstream templog;
    std::string outDir;
    std::map<std::string, int> warnCount;
};

class Loggable {
public:
    explicit Loggable(EasyLog* log = 0)
        : m_log(log)
    {
    }
    EasyLog* GetLogger() const { return m_log; }
    void SetLogger(EasyLog* log) { m_log = log; }
protected:
    EasyLog* m_log;
};

/**
 * Convert almost anything to a string
 */
template <typename type>
inline std::string stringify(type from)
{
    std::ostringstream s;
    if (!(s << from))
        throw std::logic_error("Stringify failed");
    return s.str();
}

/**
 * Output a vector<int>
 */
inline std::ostream& operator<<(std::ostream& out, std::vector<int> x)
{
    out << "[ ";
    for (unsigned i = 0; i < x.size(); i++)
        out << x[i] << " ";
    out << "]";
    return out;
}
