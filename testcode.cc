/*  skeletonapp.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2003 University of Oxford  */

/*  CCOPYRIGHT  */

// Skeleton application framework for using newimage

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"
#include "warpfns/warpfns.h"
#include "Update_deformation.h"

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="skeletonapp (Version 1.0)\nCopyright(c) 2003, University of Oxford (Mark Jenkinson)";
string examples="skeletonapp [options] --in1=<image1> --in2=<image2>";

// Each (global) object below specificies as option and can be accessed
//  anywhere in this file (since they are global).  The order of the
//  arguments needed is: name(s) of option, default value, help message,
//       whether it is compulsory, whether it requires arguments
// Note that they must also be included in the main() function or they
//  will not be active.

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<string> inname(string("-i"), string(""),
		  string("input image"),
		  true, requires_argument);
Option<string> refname(string("-r"), string(""),
		  string("reference image"),
		  true, requires_argument);
Option<string> outname(string("-o"), string(""),
		  string("reference image"),
		  true, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Local functions

// for example ... print difference of COGs between 2 images ...
int do_work(int argc, char* argv[]) 
{
  volume4D<float> vin, vref, vout, defx, defy, defz;
  volume4D<float> pdefx, pdefy, pdefz;
  read_volume4D(vin,inname.value());
  read_volume4D(vref,refname.value());
  defx=vin*0.0f;
  defy=defx;
  defz=defx;
  pdefx=defx;
  pdefy=defx;
  pdefz=defx;
  // testing apply_warp stuff
  //apply_warp(vin[0],vref[0],defx);
  UpdateDeformation(vin, vref, 2, pdefx, pdefy, pdefz, vout, defx, defy, defz);
  print_volume_info(vout);
  save_volume4D(vout,outname.value());
  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(refname);
    options.add(outname);
    options.add(verbose);
    options.add(help);
    
    nonoptarg = options.parse_command_line(argc, argv);

    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  // Call the local functions

  return do_work(argc,argv);
}

