/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

 Adrian Groves and Michael Chappell, FMRIB Image Analysis & IBME QuBIc groups

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabber_core.h"
#include "fwdmodel.h"
#include "inference.h"
#include "rundata_newimage.h"
#include "version.h"

#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <typeinfo>
#include <vector>

using namespace std;

/**
 * Print version information.
 */
static void Version()
{
    cout << "Fabber " << fabber_release_version() <<  " Source ref: " << fabber_source_version() << ", " << fabber_source_date() << endl;
}

/**
 * Print usage information.
 */
static void Usage()
{
    Version();
    cout << "Usage: fabber [--<option>|--<option>=<value> ...]" << endl
         << endl
         << "Use -f <file> to read options in option=value form" << endl
         << "Use -@ <file> to read options in command line form (DEPRECATED)." << endl
         << endl
         << "General options " << endl
         << endl;

    vector<OptionSpec> options;
    FabberRunData::GetOptions(options);

    for (unsigned int i = 0; i < options.size(); i++)
    {
        cout << options[i] << endl;
    }
}


#ifdef _WIN32
static int setenv(const char *name, const char *value, int overwrite)
{
    if (!overwrite)
    {
        size_t envsize = 0;
        int errcode = getenv_s(&envsize, NULL, 0, name);
        if (errcode || envsize)
            return errcode;
    }
    return _putenv_s(name, value);
}
#endif

static void set_environment()
{
    // This variable needs to be set, and might be missing if FSL is not installed
    setenv("FSLOUTPUTTYPE", "NIFTI_GZ", 0);
}

/**
 * Run the default command line program
 */
int execute(int argc, char **argv)
{
    EasyLog log;
    bool gzLog = false;
    int ret = 1;

    try
    {
        set_environment();

        // Create a new Fabber run
        FabberRunDataNewimage *params = new FabberRunDataNewimage(true);
        params->SetLogger(&log);
        params->Parse(argc, argv);

        // Print usage information if no arguments given, or
        // if --help specified
        if (params->GetBool("help") || argc == 1)
        {
            string model = params->GetStringDefault("model", "");
            string method = params->GetStringDefault("method", "");
            if (model != "")
            {
                FwdModel::UsageFromName(model, cout);
            }
            else if (method != "")
            {
                InferenceTechnique::UsageFromName(method, cout);
            }
            else
            {
                Usage();
            }

            return 0;
        }
        if (params->GetBool("version"))
        {
            Version();
            return 0;
        }
        else if (params->GetBool("listmodels"))
        {
            vector<string> models = FwdModel::GetKnown();
            vector<string>::iterator iter;
            for (iter = models.begin(); iter != models.end(); ++iter)
            {
                cout << *iter << endl;
            }

            return 0;
        }
        else if (params->GetBool("listmethods"))
        {
            vector<string> infers = InferenceTechnique::GetKnown();
            vector<string>::iterator iter;
            for (iter = infers.begin(); iter != infers.end(); ++iter)
            {
                cout << *iter << endl;
            }

            return 0;
        }

        // Make sure command line tool creates a parameter names file
        params->SetBool("dump-param-names");

        cout << "----------------------" << endl;
        cout << "Welcome to FABBER v" << fabber_release_version() << endl;
        cout << "----------------------" << endl;
        cout << "Source ref: " << fabber_source_version() << ", " << fabber_source_date() << endl;

        log.StartLog(params->GetOutputDir());
        cout << "Logfile started: " << log.GetOutputDirectory() << "/logfile" << endl;

        params->SetExtentFromData();
        PercentProgressCheck percent;
        params->Run(&percent);

        log.ReissueWarnings();

        // Only Gzip the logfile if we exit normally
        gzLog = params->GetBool("gzip-log");
        ret = 0;
    }
    catch (const exception &e)
    {
        log.ReissueWarnings();
        log.LogStream() << "Exception caught in fabber:\n  " << e.what() << endl;
        cerr << "Exception caught in fabber:\n  " << e.what() << endl;
    }
    catch (NEWMAT::Exception &e)
    {
        log.ReissueWarnings();
        log.LogStream() << "NEWMAT exception caught in fabber:\n  " << e.what() << endl;
        cerr << "NEWMAT exception caught in fabber:\n  " << e.what() << endl;
    }
    catch (...)
    {
        log.ReissueWarnings();
        log.LogStream() << "Some other exception caught in fabber!" << endl;
        cerr << "Some other exception caught in fabber!" << endl;
    }

    if (log.LogStarted())
    {
        cout << endl
             << "Final logfile: " << log.GetOutputDirectory() << (gzLog ? "/logfile.gz" : "/logfile") << endl;
        log.StopLog(gzLog);
    }
    else
    {
        // Flush any errors to stdout as we didn't get as far as starting the logfile
        log.StartLog(cout);
        log.StopLog();
    }
    return ret;
}
