/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

 Adrian Groves and Michael Chappell, FMRIB Image Analysis & IBME QuBIc groups

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabber_core.h"
#include "fwdmodel.h"
#include "inference.h"
#include "rundata_newimage.h"
#include "version.h"
#include "tools.h"

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
    cout << "Fabber " << fabber_version() << " : " << fabber_source_date() << endl;
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
    bool simple_output = false;
    int ret = 1;

    try
    {
        set_environment();

        // Create a new Fabber run
        FabberRunDataNewimage paramso(true);
        FabberRunDataNewimage *params = &paramso;
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
            string model_name = params->GetStringDefault("model", "");
            if (model_name != "")
            {
                std::auto_ptr<FwdModel> model(FwdModel::NewFromName(model_name));
                cout << model->ModelVersion() << endl;
            }
            else
            {
                Version();
            }
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
        else if (params->GetBool("listparams")) 
        {
            string model = params->GetStringDefault("model", "");
            std::auto_ptr<FwdModel> fwd_model(FwdModel::NewFromName(model));
            EasyLog log;
            fwd_model->SetLogger(&log); // We ignore the log but this stops it going to cerr
            fwd_model->Initialize(*params);

            vector<Parameter> model_params;
            fwd_model->GetParameters(*params, model_params);
            vector<Parameter>::iterator iter;
            for (iter = model_params.begin(); iter != model_params.end(); ++iter)
            {
                cout << iter->name << endl;
            }

            return 0;
        }
        else if (params->GetBool("listoutputs")) 
        {
            string model = params->GetStringDefault("model", "");
            std::auto_ptr<FwdModel> fwd_model(FwdModel::NewFromName(model));
            EasyLog log;
            fwd_model->SetLogger(&log); // We ignore the log but this stops it going to cerr
            fwd_model->Initialize(*params);

            vector<string> model_outputs;
            fwd_model->GetOutputs(model_outputs);
            vector<string>::iterator iter;
            for (iter = model_outputs.begin(); iter != model_outputs.end(); ++iter)
            {
                cout << *iter << endl;
            }

            return 0;
        }
        else if (params->GetBool("evaluate")) 
        {
            string model = params->GetStringDefault("model", "");
            std::auto_ptr<FwdModel> fwd_model(FwdModel::NewFromName(model));
            EasyLog log;
            fwd_model->SetLogger(&log); // We ignore the log but this stops it going to cerr
            fwd_model->Initialize(*params);

            Matrix param_values = fabber::read_matrix_file(params->GetString("evaluate-params"));
            ColumnVector p_vec = param_values.Column(1);

            int n_ts = params->GetInt("evaluate-nt", 0);
            ColumnVector data_vec(n_ts);
            if (params->HaveKey("evaluate-data")) {
                Matrix data_values = fabber::read_matrix_file(params->GetString("evaluate-data"));
                data_vec = data_values.Column(0);
            }
            else {
                data_vec = 0;
            }

            ColumnVector coords(3);
            coords(1) = 1;
            coords(2) = 1;
            coords(3) = 1;
            fwd_model->PassData(1, data_vec, coords);

            ColumnVector o_vec(n_ts);
            fwd_model->EvaluateModel(p_vec, o_vec, params->GetStringDefault("evaluate", ""));
            for (unsigned int i = 0; i < o_vec.Nrows(); i++)
            {
                cout << o_vec(i+1) << endl;
            }

            return 0;
        }
        // Make sure command line tool creates a parameter names file
        params->SetBool("dump-param-names");
        // Link to latest run
        params->SetBool("link-to-latest");
        params->SetExtentFromData();
        simple_output = params->GetBool("simple-output");

        log.StartLog(params->GetOutputDir());
        if (!simple_output)
        {
            cout << "----------------------" << endl;
            cout << "Welcome to FABBER " << fabber_version() << endl;
            cout << "----------------------" << endl;
            cout << "Last commit: " << fabber_source_date() << endl;
            cout << "Logfile started: " << log.GetOutputDirectory() << "/logfile" << endl;
            PercentProgressCheck progress;
            params->Run(&progress);
        }
        else
        {
            SimpleProgressCheck progress;
            params->Run(&progress);
        }

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
        if (!simple_output)
        {
            cout << endl
                << "Final logfile: " << log.GetOutputDirectory()
                << (gzLog ? "/logfile.gz" : "/logfile") << endl;
        }
        log.StopLog(gzLog);
    }
    else
    {
        // Flush any errors to stderr as we didn't get as far as starting the logfile
        log.StartLog(cerr);
        log.StopLog();
    }
    return ret;
}
