/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include <iostream>
#include <exception>
#include <stdexcept>
#include <map>
#include <string>
#include "inference.h"
using namespace std;
using namespace MISCMATHS;

#include "easylog.h"
using namespace Utilities;

/*** Function declarations ***/

void Usage(const string& errorString = "");

/*** Function implementations ***/


int main(int argc, char** argv)
{
  bool gzLog = false;
  try
    {
      cout << "------------------\n";
      cout << "Welcome to FABBER v2.0" << endl;
      //cout << "Welcome to FABBER development version (1.9)" << endl;

      EasyOptions args(argc, argv);

      if (args.ReadBool("help") || argc==1) 
        { 
            string model = args.ReadWithDefault("model","");
            if (model == "")
                Usage();
            else
                FwdModel::ModelUsageFromName(model, args);
                                 
            return 0; 
        }

      if (args.ReadBool("params"))
	{ 
	  string outputDir = args.ReadWithDefault("output",".");
	  EasyLog::StartLog(outputDir,false); 
	  ofstream paramFile(( EasyLog::GetOutputDirectory() + "/paramnames.txt").c_str());
	  vector<string> paramNames;
	  FwdModel* model;
	  model = FwdModel::NewFromName(args.Read("model"),args);
	  model->NameParams(paramNames);
	  for (unsigned i = 0; i < paramNames.size(); i++)
	    {
	      LOG << "      " << paramNames[i] << endl;
	      paramFile << paramNames[i] << endl;
	    }
	  paramFile.close();

	  return 0;
	}

      EasyLog::StartLog(
        args.Read("output", "Must specify an output directory, for example: --output=mytestrun"),
        args.ReadBool("overwrite"));

        
      LOG_ERR("Logfile started: " << EasyLog::GetOutputDirectory() 
	          << "/logfile" << endl);

      time_t startTime;
      time(&startTime);
      LOG_ERR("Start time: " << ctime(&startTime));

      // Diagnostic information: software versions
      // This only versions this file... should really use all.
//      LOG_ERR("FABBER development revision: $Id: fabber.cc,v 1.29 2013/01/25 15:34:35 chappell Exp $\n");
      LOG_ERR("FABBER release v2.0 \n");
      LOG << "Command line and effective options:\n" << args.Read("") << endl;
      LOG << "--output='" << EasyLog::GetOutputDirectory() << "'" << endl;
      LOG << args << "--------------------" << endl;

      // Start timing/tracing if requested
      bool recordTimings = false;
  
      if (args.ReadBool("debug-timings")) 
        { recordTimings = true; Tracer_Plus::settimingon(); }
      if (args.ReadBool("debug-instant-stack")) 
        { Tracer_Plus::setinstantstackon(); } // instant stack isn't used?
      if (args.ReadBool("debug-running-stack")) 
        { Tracer_Plus::setrunningstackon(); }
      gzLog = args.ReadBool("gzip-log");

      Tracer_Plus tr("FABBER main (outer)");
      // can't start it before this or it segfaults if an exception is thown with --debug-timings on.

      // Start a new tracer for timing purposes
      { Tracer_Plus tr2("FABBER main()");

      InferenceTechnique* infer = 
        InferenceTechnique::NewFromName(args.Read("method"));

      infer->Setup(args);
      infer->SetOutputFilenames(EasyLog::GetOutputDirectory());
      
      DataSet allData;
      allData.LoadData(args);

      // Arguments should all have been used by now, so complain if there's anything left.
      args.CheckEmpty();   
      
      // Calculations
      infer->DoCalculations(allData);
      infer->SaveResults(allData);
      delete infer;
      
      LOG_ERR("FABBER is all done." << endl);

      time_t endTime;
      time(&endTime);
      LOG << "Start time: " << ctime(&startTime);   // Bizarrely, ctime() ends with a \n.
      LOG << "End time: " << ctime(&endTime);
      LOG_ERR("Duration: " << endTime-startTime << " seconds." << endl);

      } // End of timings
     
      if (recordTimings) {
        tr.dump_times(EasyLog::GetOutputDirectory());
        LOG_ERR("Timing profile information recorded to " 
		<< EasyLog::GetOutputDirectory() << "/timings.html" << endl);
      }

      Warning::ReissueAll();

      cout << "Logfile was: " << EasyLog::GetOutputDirectory() << (gzLog ? "/logfile.gz" : "/logfile") << endl;
      EasyLog::StopLog(gzLog);

      return 0;
    }
  catch (const Invalid_option& e)
    {
      Warning::ReissueAll();
      LOG_ERR_SAFE("Invalid_option exception caught in fabber:\n  " << Exception::what() << endl);
      Usage(Exception::what());
    }
  catch (const exception& e)
    {
      Warning::ReissueAll();
      LOG_ERR_SAFE("STL exception caught in fabber:\n  " << e.what() << endl);
    }
  catch (Exception)
    {
      Warning::ReissueAll();
      LOG_ERR_SAFE("NEWMAT exception caught in fabber:\n  " 
	      << Exception::what() << endl);
    }
  catch (...)
    {
      Warning::ReissueAll();
      LOG_ERR_SAFE("Some other exception caught in fabber!" << endl);
    }
  
  if (EasyLog::LogStarted())
    {
      // Only gzip the logfile if we exited normally.
      cout << "Logfile was: " << EasyLog::GetOutputDirectory() << "/logfile" << endl;
      EasyLog::StopLog();

      //cout << "Logfile was: " << EasyLog::GetOutputDirectory() << (gzLog ? "/logfile.gz" : "/logfile") << endl;
      //EasyLog::StopLog(gzLog);
    }

  return 1;
}


void Usage(const string& errorString)
{
    cout << "\n\nUsage: fabber <arguments>\n"
     << "Arguments are mandatory unless they appear in [brackets].\n"
     << "Use -@ argfile to read additional arguments from a text file.\n\n";

    cout << "  [--help] : print this usage message\n"
     << "  --output=/path/to/output : put output here (including logfile)\n"
     << "  --method={vb|spatialvb} : use VB (or VB with spatial priors)\n"
     << "  [--max-iterations=NN] : number of iterations of VB to use (default: 10)\n"
     << "  [--data-order={interleave|concatenate|singlefile}] : should time points from multiple data "
     << "be interleaved (e.g. TE1/TE2) or left in order? (default: interleave)\n"
     << "  --data1=file1, [--data2=file2]. (use --data=file instead if --data-order=singlefile)\n"
     << "  --mask=maskfile : inference will only be performed where mask value > 0\n"
     << "  --model={quipss2|q2tips-dualecho|pcasl-dualecho} : forward model to use. "
     << "For model parameters use fabber --help --model=<model_of_interest>\n"
     << "  --noise={ar1|white} : Noise model to use\n"
     << "    ar1: two AR(1) models (optional cross-linking between TE1 & TE2)\n"
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

























