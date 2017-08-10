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

using namespace std;

EasyLog::EasyLog()
    : m_stream(0)
    , m_outdir("")
{
}

EasyLog::~EasyLog()
{
}
void EasyLog::StartLog(const string &outDir)
{
    assert(m_stream == NULL);
    assert(outDir != "");

    m_stream = new ofstream((outDir + "/logfile").c_str());
    if (!m_stream->good())
    {
        delete m_stream;
        m_stream = NULL;
        cerr << "Cannot open logfile in " << outDir << endl;
        throw runtime_error("Cannot open logfile!");
    }

    // Flush any temporary logging
    *m_stream << m_templog.str() << flush;
    m_outdir = outDir;
}

void EasyLog::StartLog(ostream &s)
{
    assert(m_stream == NULL);
    m_stream = &s;
    m_outdir = "";

    // Flush any temporary logging
    *m_stream << m_templog.str() << flush;
}

void EasyLog::StopLog(bool gzip)
{
    assert(m_stream != NULL);

    if (m_outdir != "")
    {
        if (gzip)
        {
#ifdef _WIN32
            LogStream() << "EasyLog::GZIP logfile not supported under Windows" << std::endl;
#else
            int retVal = system(("gzip " + m_outdir + "/logfile").c_str());
            if (retVal != 0)
                LogStream() << "Failed to gzip logfile.  Oh well." << std::endl;
#endif
        }
        // We created this ofstream and need to tidy it up
        delete m_stream;
    }

    // Leave variables in consistent state in case
    // somebody calls StartLog() again
    m_stream = NULL;
    m_outdir = "";
}

bool EasyLog::LogStarted()
{
    return m_stream != NULL;
}
const string &EasyLog::GetOutputDirectory()
{
    return m_outdir;
}
std::ostream &EasyLog::LogStream()
{
    if (m_stream == NULL)
    {
        return m_templog;
    }
    else
    {
        return *m_stream;
    }
}

void EasyLog::WarnOnce(const string &text)
{
    if (++m_warncount[text] == 1)
        LogStream() << "WARNING ONCE: " << text << std::endl;
}

void EasyLog::WarnAlways(const string &text)
{
    ++m_warncount[text];
    LogStream() << "WARNING ALWAYS: " << text << std::endl;
}

void EasyLog::ReissueWarnings()
{
    if (m_warncount.size() == 0)
        return; // avoid issuing pointless message

    LogStream() << "\nSummary of warnings (" << m_warncount.size() << " distinct warnings)\n";
    for (map<string, int>::iterator it = m_warncount.begin(); it != m_warncount.end(); ++it)
        LogStream() << "Issued "
                    << ((it->second == 1) ? "once: " : stringify(it->second) + " times: ")
                    << it->first << std::endl;
}
