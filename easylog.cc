/*  easylog.cc - a fairly minimal logging-to-file implementation

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "easylog.h"

#include "assert.h"
#include <stdexcept>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <stdlib.h>

#ifdef _WIN32
#include "direct.h"
#else
#include <sys/stat.h>
#endif

using namespace std;

static bool is_dir(string path)
{
#ifdef _WIN32
	struct _stat s;
	if (_stat(path.c_str(), &s) == 0)
#else
	struct stat s;
	if (stat(path.c_str(), &s) == 0)
#endif
	{
		return (s.st_mode & S_IFDIR);
	}
	else
	{
		// Does not exist, so not a directory...
		return false;
	}
}
#ifdef _WIN32
	map<DWORD, EasyLog> EasyLog::s_logs;
#else
  #ifdef USE_PTHREADS
	map<pthread_t, EasyLog> EasyLog::s_logs;
  #else
	EasyLog *EasyLog::s_log=0;
  #endif
#endif

EasyLog& EasyLog::CurrentLog()
{
#ifdef _WIN32
	if (s_logs.count(GetCurrentThreadId()) == 0) {
		s_logs[GetCurrentThreadId()] = EasyLog();
	}
	return s_logs[GetCurrentThreadId()];
#else
  #ifdef USE_PTHREADS
	if (s_logs.count(pthread_self()) == 0) {
		s_logs[pthread_self()] = EasyLog();
	}
	return s_logs[pthread_self()];
 #else
	if (s_log) return *s_log;
	else {
		s_log = new EasyLog();
		return *s_log;
	}
  #endif
#endif
}

EasyLog::EasyLog() :
		filestream(0), outDir("")
{
}

EasyLog::~EasyLog()
{
	// FIXME hack and not threadsafe
	if (s_log) {
		delete s_log;
		s_log = NULL;
	}
}

void EasyLog::StartLog(const string& basename, bool overwrite, bool link_to_latest)
{
	assert(filestream == NULL);
	assert(basename != "");
	outDir = basename;

	// From Wooly's utils/log.cc
	int count = 0;
	while (true)
	{
		if (count >= 50) // I'm using a lot for some things
		{
			throw std::runtime_error(
					("Cannot create directory (bad path, or too many + signs?):\n    " + outDir).c_str());
		}

		// Clear errno so it can be inspected later; result is only meaningful if mkdir fails.
		errno = 0;
		int ret = 0;

#ifdef _WIN32
		ret = _mkdir(outDir.c_str());
#else
		ret = mkdir(outDir.c_str(), 0777);
#endif

		if (ret == 0) // Success, directory created
			break;
		else if (overwrite)
		{
			if ((errno == EEXIST) && is_dir(outDir))
			{
				// If directory already exists -- that's fine.
				break;
			}
			else
				// Other error -- might be a problem!
				throw std::runtime_error(
						("Unexpected problem creating directory in --overwrite mode:\n    " + outDir).c_str());
		}

		outDir += "+";
		count++;
	}

	filestream = new ofstream((outDir + "/logfile").c_str());

	if (!filestream->good())
	{
		delete filestream;
		filestream = NULL;
		cout << "Cannot open logfile in " << outDir << endl;
		throw runtime_error("Cannot open logfile!");
	}

#ifdef _WIN32
#else
	// Might be useful for jobs running on the queue:
	system(("uname -a > " + outDir + "/uname.txt").c_str());

	if (link_to_latest)
	{
		// try to make a link to the latest version
		// If this fails, it doesn't really matter.
		system(("ln -sfn '" + outDir + "' '" + basename + "_latest'").c_str());
	}
#endif

	// Flush any temporary logging
	*filestream << templog.str() << flush;
}

void EasyLog::StartLog(ostream& s)
{
	assert(filestream == NULL);
	filestream = &s;
	outDir = "";

	// Flush any temporary logging
	*filestream << templog.str() << flush;
}

const std::string& EasyLog::GetOutputDirectory()
{
	assert(filestream != NULL);
	return outDir;
}

void EasyLog::StopLog(bool gzip)
{
	assert(filestream != NULL);

	if (outDir != "")
	{
		if (gzip)
		{
#ifdef _WIN32
			LogStream() << "EasyLog::GZIP logfile not supported under Windows" << std::endl;
#else
			int retVal = system(("gzip " + outDir + "/logfile").c_str());
			if (retVal != 0)
				LogStream() << "Failed to gzip logfile.  Oh well." << std::endl;
#endif
		}
		// We created this ofstream and need to tidy it up
		delete filestream;
	}

	// Leave variables in consistent state in case
	// somebody calls StartLog() again
	filestream = NULL;
	outDir = "";
}

bool EasyLog::LogStarted()
{
	return filestream != NULL;
}

std::ostream& EasyLog::LogStream()
{
	if (filestream == NULL)
	{
		return templog;
	}
	else
	{
		return *filestream;
	}
}

void EasyLog::WarnOnce(const string& text)
{
	if (++warnCount[text] == 1)
		LogStream() << "WARNING ONCE: " << text << std::endl;
}

void EasyLog::WarnAlways(const string& text)
{
	++warnCount[text];
	LogStream() << "WARNING ALWAYS: " << text << std::endl;
}

void EasyLog::ReissueWarnings()
{
	if (warnCount.size() == 0)
		return; // avoid issuing pointless message

	LogStream() << "\nSummary of warnings (" << warnCount.size() << " distinct warnings)\n";
	for (map<string, int>::iterator it = warnCount.begin(); it != warnCount.end(); it++)
		LogStream() << "Issued " << ( (it->second==1)? "once: " : stringify(it->second)+" times: " ) << it->first << std::endl;
}
