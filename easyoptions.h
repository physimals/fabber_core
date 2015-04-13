/*  easyoptions.h - FABBER's options-handling class

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once
#include <map>
#include <string>
#include "easylog.h"

class EasyOptions {
 public:
    // Create an instance from command-line options
    EasyOptions(int argc, char** argv);
        // break down options into key=value clauses

    EasyOptions(const map<string,string>& src, 
                const map<string,const Matrix*>* dataIn, 
                map<string,Matrix>* dataOut) 
        : args(src) // copy key=value pairs from src; use key="" for no-argument options
        { 
          assert(inMatrices==NULL); 
          assert(outMatrices==NULL); 
          inMatrices=dataIn; 
          outMatrices=dataOut; 
          assert(outMatrices->empty());
	  args[""] = "fabber_library"; // This would normally hold argv[0].
        }

    // Below: option-reading values.  Once they are called, 
    // the corresponding key=value pair is removed... this is deliberate, to
    // ensure that eacy argument is used exactly once.  If you want to use an
    // argument in more than one place, you're probably doing something wrong...
    // this should make code more maintainable by making it easy to see where
    // options have an effect on code.
    
    string Read(const string& key)
        { return Read(key, "Missing mandatory option: --" + key + "\n"); }
        // Get a mandatory option.
        // throws if option is missing or has no =value clause
    string Read(const string& key, const string& msg);
        // supply msg if you want a customized (helpful) error message.
        
    bool ReadBool(const string& key);
        // true if present, false if absent, throws error if --key=value.
        
    string ReadWithDefault(const string& key, const string& def);
        // Like default, but if --key=value is omitted it just treats
        // as if --key=default.
        
    void CheckEmpty();
        // throws if there are any options left

    ~EasyOptions() { }; 
        // throwing an exception in a destructor would be a bad idea!
        // Also, note that args may not be empty, if an exception has
        // been thrown.

    friend ostream& operator<<(ostream& out, const EasyOptions& opts);

private:
    map<string,string> args;
        // all remaining unused options.
    void AddKey(const string& key); // internal helper function

    // These are only used in fabber_library mode:
    static const map<string,const Matrix*>* inMatrices;
    static map<string,Matrix>* outMatrices;
public:
    static bool UsingMatrixIO() { assert(!inMatrices == !outMatrices); return (inMatrices!=NULL); }
    static const Matrix& InMatrix(const string& filename);  // simpler syntax
    static Matrix& OutMatrix(const string& filename);
};
 
// Helper function:
Matrix read_vest_fabber(const string& filename);

typedef EasyOptions ArgsType;

// Use NEWMAT's exception handling base class
class Invalid_option : public Runtime_error
{
public:
  static unsigned long Select;          // for identifying exception
  Invalid_option(const char* c) : Runtime_error(c) { return; }
  Invalid_option(const string& s) : Runtime_error(s.c_str()) { return; }
};

// Convert a string into almost anything, using operator>>. See the C++ FAQ-Lite: 
// http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.3
template<typename T>
inline T convertTo(const std::string& s)
{
  T x;
  istringstream i(s);
  char c;
  if (!(i>>x) || (i.get(c)))
    throw Invalid_option(
        "convertTo failed... couldn't convert string '" + s + "' into a value.\n");
  return x;
}


