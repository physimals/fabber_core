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

#include "fwdmodel.h"
#include "inference.h"
#include "fabber_core.h"

using namespace Utilities;

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
	copy(str_vector.begin(), str_vector.end(), ostream_iterator<string>(str_stream, separator));
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

	cout << "\n\nUsage: fabber <arguments>" << endl << "Use -@ argfile to read additional arguments from a text file."

	<< endl << endl << " General options " << endl << endl
			<< "  --help                     Print this usage message. When --method or --model" << endl
			<< "                             are specified in addition, display relevant model" << endl
			<< "                             / inference method options" << endl
			<< "  --output=<dir>             Directory for output files (including logfile)" << endl
			<< "  --method=<method>          Use this inference method " << endl
			<< "                             (known inference methods: " << methods << ")" << endl
			<< "  --model=<model name>       Forward model to use. Known forward models are: " << endl
			<< "                            " << fwdmodels << endl
			<< "  --data=<file>              Specify a single input data file" << endl
			<< "  --data1=<file>, --data2=<file2>" << endl
			<< "                             Specify multiple input data files (see data-order)" << endl
			<< "  --data-order=[interleave|concatenate|singlefile]  " << endl
			<< "                             If multiple data files are specified, how they will" << endl
			<< "                             be handled: concatenate = one after the other, " << endl
			<< "                             interleave = first record from each file, then " << endl
			<< "                             second, etc." << endl
			<< "  --mask=<file>              Mask file. Inference will only be performed where " << endl
			<< "                             mask value > 0" << endl
			<< "  --save-model-fit           Save a file containing the mean model prediction " << endl
			<< "                             at each voxel" << endl << "  --save-residuals           " << endl
			<< "  --allow-bad-voxels         Skip to next voxel if a numerical exception occurs" << endl << endl
			<< " Noise options " << endl << endl << "  --ar1-cross-terms={dual|same|none}]  " << endl
			<< "                             Two types of cross-linking, or none (default: dual)" << endl
			<< "  --noise-pattern=<phi_index_pattern>] " << endl
			<< "                             repeating pattern of noise variances for each point" << endl
			<< "                             (e.g. --noise-pattern=12 gives odd and even data " << endl
			<< "                             points different noise variances)" << endl << endl << " VB options "
			<< endl << endl << "  --noise=<noise model name> Noise model to use. Known noise models are: " << endl
			<< "                            " << noisemodels << endl
			<< "  --max-iterations=<n>       number of iterations of VB to use (default: 10)" << endl
			<< "  --print-free-energy        Calculate & dump F to the logfile after each update" << endl << endl
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
	if (errorString.length() > 0)
		cout << "\nImmediate cause of error: " << errorString << endl;
}

