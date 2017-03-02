/*  easylog.cc - a fairly minimal logging-to-file implementation

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "easylog.h"

#include <assert.h>
#include <errno.h>
#include <fstream>
#include <stdexcept>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

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
    } else {
        // Does not exist, so not a directory...
        return false;
    }
}

EasyLog::EasyLog()
    : stream(0)
    , outDir("")
{
}

EasyLog::~EasyLog()
{
}

void EasyLog::StartLog(const string& outDir)
{
    assert(stream == NULL);
    assert(outDir != "");

    stream = new ofstream((outDir + "/logfile").c_str());
    if (!stream->good()) {
        delete stream;
        stream = NULL;
        cerr << "Cannot open logfile in " << outDir << endl;
        throw runtime_error("Cannot open logfile!");
    }

    // Flush any temporary logging
    *stream << templog.str() << flush;
}

void EasyLog::StartLog(ostream& s)
{
    assert(stream == NULL);
    stream = &s;
    outDir = "";

    // Flush any temporary logging
    *stream << templog.str() << flush;
}

void EasyLog::StopLog(bool gzip)
{
    assert(stream != NULL);

    if (outDir != "") {
        if (gzip) {
#ifdef _WIN32
            LogStream() << "EasyLog::GZIP logfile not supported under Windows" << std::endl;
#else
            int retVal = system(("gzip " + outDir + "/logfile").c_str());
            if (retVal != 0)
                LogStream() << "Failed to gzip logfile.  Oh well." << std::endl;
#endif
        }
        // We created this ofstream and need to tidy it up
        delete stream;
    }

    // Leave variables in consistent state in case
    // somebody calls StartLog() again
    stream = NULL;
    outDir = "";
}

bool EasyLog::LogStarted()
{
    return stream != NULL;
}

const string& EasyLog::GetOutputDirectory()
{
    return outDir;
}

std::ostream& EasyLog::LogStream()
{
    if (stream == NULL) {
        return templog;
    } else {
        return *stream;
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
    for (map<string, int>::iterator it = warnCount.begin(); it != warnCount.end(); ++it)
        LogStream() << "Issued " << ((it->second == 1) ? "once: " : stringify(it->second) + " times: ") << it->first << std::endl;
}
