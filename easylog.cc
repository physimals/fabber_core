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

ostream* EasyLog::filestream = NULL;
string EasyLog::outDir = "";
stringstream EasyLog::templog;

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
			LOG << "want to overwrite" << endl;
			if ((errno == EEXIST) && is_dir(outDir))
			{
				// If directory already exists -- that's fine.
				LOG << "already there, fine" << endl;
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

void EasyLog::StopLog(bool gzip)
{
	assert(filestream != NULL);

	if (outDir != "")
	{
		if (gzip)
		{
#ifdef _WIN32
			LOG << "EasyLog::GZIP logfile not supported under Windows" << std::endl;
#else
			int retVal = system(("gzip " + outDir + "/logfile").c_str());
			if (retVal != 0)
				LOG << "Failed to gzip logfile.  Oh well." << std::endl;
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

map<string, int> Warning::issueCount;

void Warning::IssueOnce(const string& text)
{
	if (++issueCount[text] == 1)
		LOG_ERR("WARNING ONCE: " << text << std::endl);
}

void Warning::IssueAlways(const string& text)
{
	++issueCount[text];
	LOG_ERR("WARNING ALWAYS: " << text << std::endl);
}

void Warning::ReissueAll()
{
	if (issueCount.size() == 0)
		return; // avoid issuing pointless message

	LOG_ERR("\nSummary of warnings (" << issueCount.size() << " distinct warnings)\n");
	for (map<string, int>::iterator it = issueCount.begin(); it != issueCount.end(); it++)
		LOG_ERR(
				"Issued " << ( (it->second==1)? "once: " : stringify(it->second)+" times: " ) << it->first << std::endl);
}
