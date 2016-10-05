/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

 Adrian Groves and Michael Chappell, FMRIB Image Analysis & IBME QuBIc groups

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabber_core.h"

#include "fwdmodel.h"
#include "inference.h"

#include "utils/tracer_plus.h"

#include <exception>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <memory>
#include <vector>

using Utilities::Tracer_Plus;

static OptionSpec OPTIONS[] =
		{
				{ "help", OPT_BOOL,
						"Print this usage method. If given with --method or --model, display relevant method/model usage information",
						OPT_NONREQ, "" },
				{ "listmethods", OPT_BOOL, "List all known inference methods", OPT_NONREQ, "" },
				{ "listmodels", OPT_BOOL, "List all known forward models", OPT_NONREQ, "" },
				{ "output", OPT_STR, "Directory for output files (including logfile)", OPT_REQ, "" },
				{ "overwrite", OPT_BOOL,
						"If set will overwrite existing output. If not set, new output directories will be created by appending '+' to the directory name ",
						OPT_NONREQ, "" },
				{ "link-to-latest", OPT_BOOL,
						"If set will try to create a link to the most recent output directory with the prefix _latest",
						OPT_NONREQ, "" },
				{ "method", OPT_STR, "Use this inference method", OPT_NONREQ, "" },
				{ "model", OPT_STR, "Use this forward model", OPT_NONREQ, "" },
				{ "data", OPT_FILE, "Specify a single input data file", OPT_REQ, "" },
				{ "data<n>", OPT_FILE, "Specify multiple data files for n=1, 2, 3...", OPT_NONREQ, "" },
				{ "data-order", OPT_STR,
						"If multiple data files are specified, how they will be handled: concatenate = one after the other,  interleave = first record from each file, then  second, etc.",
						OPT_NONREQ, "interleave" },
				{ "mask", OPT_FILE, "Mask file. Inference will only be performed where mask value > 0", OPT_NONREQ, "" },
				{ "save-model-fit", OPT_BOOL, "Save the model prediction as a 4d volume", OPT_NONREQ, "" },
				{ "save-residuals", OPT_BOOL,
						"Save the difference between the data and the model prediction as a 4d volume", OPT_NONREQ, "" },
				{ "" }, };

/**
 * Print usage information.
 */
void Usage()
{
	cout << "\n\nUsage: fabber [--<option>|--<option>=<value> ...]" << endl << endl
			<< "Use -@ <file> to read additional arguments in command line form from a text file (DEPRECATED)." << endl
			<< "Use -f <file> to read options in option=value form" << endl << endl << "General options " << endl
			<< endl;

	for (int i = 0; OPTIONS[i].name != ""; i++)
	{
		cout << OPTIONS[i] << endl;
	}
#if 0
	<< " Spatial VB options: " << endl << endl << "  --param-spatial-priors=<choice_of_prior_forms> " << endl
	<< "                             Specify a type of prior to use for each forward " << endl
	<< "                             model parameter.  One letter per parameter.  " << endl
	<< "                             S=spatial, N=nonspatial, D=Gaussian-process-based" << endl
	<< "                             combined prior" << endl << "  --fwd-initial-prior=<prior_vest_file> "
	<< endl << "                             Specify the nonspatial prior distributions on the" << endl
	<< "                             forward model parameters.  The vest file is the " << endl
	<< "                             covariance matrix supplemented by the prior means" << endl
	<< "                             See the documentation for details.  Very important" << endl
	<< "                             if 'D' prior is used." << endl << endl;
#endif
}

/**
 * Run the default command line program
 */
int execute(int argc, char** argv)
{
	bool gzLog = false;
	try
	{
		// Create a new Fabber run
		FabberRunData params;
		PercentProgressCheck percent;
		params.SetProgressCheck(&percent);
		params.Parse(argc, argv);

		// Print usage information if no arguments given, or
		// if --help specified
		if (params.GetBool("help") || argc == 1)
		{
			string model = params.GetStringDefault("model", "");
			string method = params.GetStringDefault("method", "");
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
		else if (params.GetBool("listmodels"))
		{
			vector < string > models = FwdModel::GetKnown();
			vector<string>::iterator iter;
			for (iter = models.begin(); iter != models.end(); iter++)
			{
				cout << *iter << endl;
			}

			return 0;
		}
		else if (params.GetBool("listmethods"))
		{
			vector < string > infers = InferenceTechnique::GetKnown();
			vector<string>::iterator iter;
			for (iter = infers.begin(); iter != infers.end(); iter++)
			{
				cout << *iter << endl;
			}

			return 0;
		}

		cout << "----------------------" << endl;
		cout << "Welcome to FABBER v" << VERSION << endl;
		cout << "----------------------" << endl;

		EasyLog::StartLog(params.GetStringDefault("output", "."), params.GetBool("overwrite"), params.GetBool("link-to-latest"));
		cout << "Logfile started: " << EasyLog::GetOutputDirectory() << "/logfile" << endl;

		// Start timing/tracing if requested
		bool recordTimings = false;

		if (params.GetBool("debug-timings"))
		{
			recordTimings = true;
			Tracer_Plus::settimingon();
		}
		if (params.GetBool("debug-instant-stack"))
		{
			Tracer_Plus::setinstantstackon();
		} // instant stack isn't used?
		if (params.GetBool("debug-running-stack"))
		{
			Tracer_Plus::setrunningstackon();
		}

		// can't start it before this or it segfaults if an exception is thown with --debug-timings on.
		Tracer_Plus tr("FABBER main (outer)");

		// Start a new tracer for timing purposes
		{
			params.Run();
		}

		if (recordTimings)
		{
			tr.dump_times(EasyLog::GetOutputDirectory());
			LOG << "Timing profile information recorded to " << EasyLog::GetOutputDirectory() << "/timings.html"
					<< endl;
		}

		Warning::ReissueAll();

		gzLog = params.GetBool("gzip-log");
		cout << "Final logfile: " << EasyLog::GetOutputDirectory() << (gzLog ? "/logfile.gz" : "/logfile") << endl;
		EasyLog::StopLog(gzLog);

		return 0;
	} catch (const DataNotFound& e)
	{
		Warning::ReissueAll();
		LOG_ERR("Data not found:\n  " << e.what() << endl);
	} catch (const Invalid_option& e)
	{
		Warning::ReissueAll();
		LOG_ERR("Invalid_option exception caught in fabber:\n  " << e.what() << endl);
	} catch (const exception& e)
	{
		Warning::ReissueAll();
		LOG_ERR("STL exception caught in fabber:\n  " << e.what() << endl);
	} catch (Exception)
	{
		Warning::ReissueAll();
		LOG_ERR("NEWMAT exception caught in fabber:\n  " << Exception::what() << endl);
	} catch (...)
	{
		Warning::ReissueAll();
		LOG_ERR("Some other exception caught in fabber!" << endl);
	}

	if (EasyLog::LogStarted())
	{
		// FIXME Gzip the logfile if we exited normally?
		cout << "Logfile was: " << EasyLog::GetOutputDirectory() << "/logfile" << endl;
		EasyLog::StopLog();
	}
	else
	{
		// Flush any errors to stdout as we didn't get as far as starting the logfile
		EasyLog::StartLog (cout);
	}
	return 1;
}
