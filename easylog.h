/*  easylog.h - a fairly minimal logging-to-file implementation

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "utils/tracer_plus.h"
#include <iostream>
#include <string>
#include "assert.h"

using namespace std;
using namespace Utilities;

#define PRINTNOTE fprintf(stderr, "Note: %s line %d\n", __FILE__, __LINE__);

#define LOG (*EasyLog::CurrentLog())
// Use LOG just like you'd use cout.  (i.e. LOG << your_data << endl;)
// Make sure you call StartLog or StartLogUsingStream before LOGging anything.
// If you want to dump variables directly to this stream, just implement
//    ostream& operator<<(ostream& stream, const your_type& data);

// Two new ones that aren't used very much
#define LOG_SAFE_ELSE_CERR(x) (EasyLog::LogStarted()?(void)(LOG<<x):(void)(cout<<x)) 
#define LOG_SAFE_ELSE_DISCARD(x) (EasyLog::LogStarted()?(void)(LOG<<x):(void)(0)) 


#ifdef __FABBER_LIBRARYONLY
// These names no longer have their original meaning in the library version
#define LOG_ERR(x) ((void)((LOG<<x))) 
#define LOG_ERR_SAFE(x) LOG_ERR(x)
#define LOG_SAFE_ELSE_CERR(x) (EasyLog::LogStarted()?(void)(LOG<<x):(void)(cout<<x)) 
#define LOG_SAFE_ELSE_DISCARD(x) (EasyLog::LogStarted()?(void)(LOG<<x):(void)(0)) 

#else //__FABBER_LIBRARYONLY

// Changed these to use cout rather than cerr -- since the logfile contains 
// almost everything I don't expect people to redirect cout very often.

#define LOG_ERR(x) ((void)((cout<<x),(LOG << x)))
// Use only for things that are safe to duplicate
// LOG_ERR("Note: the matrix " << name << " == " << mat.t() << endl);

#define LOG_ERR_SAFE(x) (cout<<x,EasyLog::LogStarted()?(void)(LOG<<x):(void)0)
// This is safe to use if even if the log might not have started.
// It should only be used in main()'s exception-catching routines


#endif //__FABBER_LIBRARYONLY

class EasyLog {
 public:
  static ostream* CurrentLog()
    { assert(filestream != NULL); return filestream; }
  static const string& GetOutputDirectory()
    { assert(filestream != NULL); return outDir; }

  static void StartLog(const string& basename, bool overwrite);
  static void StartLogUsingStream(ostream& s);
  static void StopLog(bool gzip = false);

  static bool LogStarted() 
    { return filestream != NULL; }
  // only use this in situations where the log might not have been started..
  // e.g. in main()'s exception-handling routines

 private:
  static ostream* filestream;
  static string outDir;
};

// Other useful functions:

// Convert almost anything to a string
#include <sstream>
#include <stdexcept>
template<typename type>
inline string stringify(type from)
{
  ostringstream s;
  if (!(s << from))
    throw logic_error("Stringify failed");
  return s.str();
}



// output a vector<int>
#include <vector>
inline ostream& operator<<(ostream& out, vector<int> x)
{
  out << "[ ";
  for (unsigned i = 0; i < x.size(); i++)
    out << x[i] << " ";
  out << "]";
  return out;
}

#include <map>

// A work in progress:
class Warning {
 public:
  static void IssueOnce(const string& text);
  static void IssueAlways(const string& text);
  static void ReissueAll();
 private:
  static map<string, int> issueCount;
};
