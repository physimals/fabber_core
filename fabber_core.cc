/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

 Adrian Groves and Michael Chappell, FMRIB Image Analysis & IBME QuBIc groups

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include <exception>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <memory>
#include <vector>

#include "fabber_core.h"

void libfabber_present(void)
{
}

/**
 * Run the default command line program
 */
int execute(int argc, char** argv)
{
	bool gzLog = false;
	try
	{
		cout << "----------------------" << endl;
		cout << "Welcome to FABBER v" << VERSION << endl;
		cout << "----------------------" << endl;

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
			if (model == "")
				Usage();
			else
			{
				FwdModel::UsageFromName(model, cout);
			}
			return 0;
		}

		EasyLog::StartLog(params.GetStringDefault("output", "."), true);
		cout << "Logfile started: " << EasyLog::GetOutputDirectory() << "/logfile" << endl;
		string outputDir = params.GetStringDefault("output", ".");

		// Diagnostic information: software versions
		LOG << "FABBER release v" << VERSION << endl;

		time_t startTime;
		time(&startTime);
		LOG << "Start time: " << ctime(&startTime);

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

			LOG << "FABBER is all done." << endl;

			time_t endTime;
			time(&endTime);
			LOG << "Start time: " << ctime(&startTime); // Bizarrely, ctime() ends with a \n.
			LOG << "End time: " << ctime(&endTime);
			LOG << "Duration: " << endTime - startTime << " seconds." << endl;

		} // End of timings

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
		// Only gzip the logfile if we exited normally.
		cout << "Logfile was: " << EasyLog::GetOutputDirectory() << "/logfile" << endl;
		EasyLog::StopLog();
	}

	return 1;
}

/**
 * Concatenate a vector of strings into a single string.
 * @param str_vector Vector of strings.
 * @param separator Separator for each string in the new string.
 * @return single string.
 */
string vectorToString(vector<string> str_vector, const char* separator)
{
	stringstream str_stream;
	copy(str_vector.begin(), str_vector.end(), ostream_iterator<string> (str_stream, separator));
	string str = str_stream.str();
	// Trim trailing delimiter.
	if (str.size() > 0)
	{
		str.resize(str.size() - 1);
	}
	return str;
}

/**
 * Print usage information.
 * @param errorString Optional error string.
 */
void Usage(const string& errorString)
{
	string fwdmodels = vectorToString(FwdModelFactory::GetInstance()->GetNames(), "|");
	string methods = vectorToString(InferenceTechniqueFactory::GetInstance()->GetNames(), "|");
	string noisemodels = vectorToString(NoiseModelFactory::GetInstance()->GetNames(), "|");

	cout << "\n\nUsage: fabber <arguments>\n" << "Arguments are mandatory unless they appear in [brackets].\n"
			<< "Use -@ argfile to read additional arguments from a text file.\n\n";
	cout << "  [--help] : print this usage message\n"
			<< "  --output=/path/to/output : put output here (including logfile)\n" << "  --method={" << methods
			<< "} : use VB (or VB with spatial priors)\n"
			<< "  [--max-iterations=NN] : number of iterations of VB to use (default: 10)\n"
			<< "  [--data-order={interleave|concatenate|singlefile}] : should time points from multiple data "
			<< "be interleaved (e.g. TE1/TE2) or left in order? (default: interleave)\n"
			<< "  --data1=file1, [--data2=file2]. (use --data=file instead if --data-order=singlefile)\n"
			<< "  --mask=maskfile : inference will only be performed where mask value > 0\n" << "  --model={"
			<< fwdmodels << "} : forward model to use. "
			<< "For model parameters use fabber --help --model=<model_of_interest>\n" << "  --noise={" << noisemodels
			<< "} : Noise model to use\n" << "    ar1: two AR(1) models (optional cross-linking between TE1 & TE2)\n"
			<< "      [--ar1-cross-terms={dual|same|none}] : two types of cross-linking, or none (default: dual)\n"
			<< "    white: white noise model, optionally with different noise variances at some data points\n"
			<< "      [--noise-pattern=<phi_index_pattern>] : repeating pattern of noise variances for each data point "
			<< "(e.g. --noise-pattern=12 gives odd and even data points different noise variances)\n"
			<< "  [--save-model-fit] and [--save-residuals] : Save model fit/residuals files\n"
			<< "  [--print-free-energy] : Calculate & dump F to the logfile after each update\n"
			<< "  [--allow-bad-voxels] : Skip to next voxel if a numerical exception occurs (don't stop)\n"
			<< "For spatial priors (using --method=spatialvb):\n"
			<< "  --param-spatial-priors=<choice_of_prior_forms>: Specify a type of prior to use for each"
			<< " forward model parameter.  One letter per parameter.  S=spatial, N=nonspatial, D=Gaussian-process-based combined prior\n"
			<< "  --fwd-initial-prior=<prior_vest_file>: specify the nonspatial prior distributions on the forward model parameters.  The vest file is the covariance matrix supplemented by the prior means; see the documentation for details.  Very important if 'D' prior is used.\n"
			<< endl;
	if (errorString.length() > 0)
		cout << "\nImmediate cause of error: " << errorString << endl;
}

