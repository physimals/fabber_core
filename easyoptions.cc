/*  easyoptions.cc - FABBER's options-handling class

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include <stdexcept> 
#include "easyoptions.h"
using namespace Utilities;
#include <fstream>
#include <sstream>
#include "miscmaths/miscmaths.h"
using namespace MISCMATHS;

// Why yes, these are basically global variables.  They're only used in fabber_library mode, and they are only accessed
// by a few select file-loading/saving functions that are friends of EasyOptions.
const map<string,const Matrix*>* EasyOptions::inMatrices = NULL;
map<string,Matrix>* EasyOptions::outMatrices = NULL;

string unescapeFilename(const string& filename)
{
   if (filename[0]=='<' && filename[filename.size()-1]=='>' && filename.size()>2)
	return filename.substr(1,filename.size()-2); // Simple "<basisSet1>" format for input matrices
   if (filename.substr(0,3) == "<>/" && filename.size()>3)
	return filename.substr(3); // Hacky "<>/output1" form for output matrices
   return ""; // unambigously means that this wasn't an escaped filename; i.e. it's presumably a real filename.
}

bool isEscapedFilename(const string& filename)
{
   bool b = unescapeFilename(filename) != "";
   return b;
}


Matrix read_vest_fabber(const string& filename)
{
   Tracer_Plus("read_vest_fabber");
   if (isEscapedFilename(filename))
     {
	return EasyOptions::InMatrix(unescapeFilename(filename));
     }
   else
     {
        return read_vest(filename);
     }
}


const Matrix& EasyOptions::InMatrix(const string& filename)
{
   if (!isEscapedFilename(filename))
      throw Invalid_option("Should be an escaped matrix name: " + filename);
   string s = unescapeFilename(filename);
   if (inMatrices->find(s) == inMatrices->end())
      throw Invalid_option("Missing input matrix " + s);
   else
      return *(inMatrices->find(s)->second);
}
   
Matrix& EasyOptions::OutMatrix(const string& filename) 
{ 
   if (!isEscapedFilename(filename))
      throw Invalid_option("Should be an escaped matrix name: " + filename);
   return (*outMatrices)[unescapeFilename(filename)]; 
}

EasyOptions::EasyOptions(int argc, char** argv) 
{
    Tracer_Plus tr("EasyOptions::EasyOptions");
    // Parse argv into key-value pairs
    // Accepted forms:
    // --key=value -> args[key] == value
    // --key value -> args[key] == value
    // --key       -> args[key] == "" 
    // For the last one, use args.count(key) as a boolean flag.  Do not look at
    //   args[key] because that will create the empty key for you!
    // Also, the command name is stored in args[""].

    // Extensions: 
    // - use a multimap so duplicate options are okay
    // - accept short-form options (e.g. -k)

    // For completeness, store the command line as well
    // Easy enough to dispose of (if you remember) 
    args[""] = argv[0];

    for (int a = 1; a < argc; a++)
    {
        if (string(argv[a],0,2) == "--")
        {
            string key = argv[a] + 2; // skip the "--"
            AddKey(key);
        }
	else if (string(argv[a]) == "-@")
        {
//            cout << "Reading additional -@ options from " << argv[a+1] << endl;
            ifstream is(argv[++a]);
            if (!is.good()) {
                cout << "Couldn't read -@ input file: '" << string(argv[a]) << "'";
                throw Invalid_option("Couldn't read input file: -@ " + string(argv[a]));
            }
            char c;
            string param;
            while (is.good())
            {
                is.get(c);
                if (!isspace(c))
                    param += c;
                else if (param == "")
		    { } // repeated whitespace, so do nothing
		else if (string(param, 0, 2) == "--") // we have an option
                { 
                    AddKey(string(param,2,string::npos)); 
                    param = ""; 
                }
                else if (string(param, 0, 1) == "#") // comment
                {
                    // discard this word and the rest of the line.
                    param = "";
                    while (is.good() && c != '\n') 
                        is.get(c);
                }
                else if (string(param, 0, 2) == "-@")
                {
                    throw Invalid_option("Sorry, at the moment you can only use -@ on the command line (no recursion).\n");
                }
                else
                {
                    throw Invalid_option("Invalid option '" + param + "' found in file '" + argv[a] + "'\n");
                }
            }       
        }
        else
        {
            cout << "Option doesn't begin with --: " + string(argv[a]);
            throw Invalid_option("An option doesn't begin with --\n");
        }
    }
}

void EasyOptions::AddKey(const string& key)
{
//    cout << "Adding key '" << key << "'" << endl;
    string::size_type eqPos = key.find("=");

    if (args.count(string(key,0,eqPos)) > 0)
        throw Invalid_option("Duplicated option: '" + key + "'\n");
    else if (eqPos != (key.npos))
        args[string(key,0,eqPos)] = string(key, eqPos+1);
//    else if (a+1 < argc && string(argv[a+1],0,1) != "-")
//        args[key] = argv[++a]; // This is the --option value rather than --option=value format.
    else
        args[key] = "";
}

string EasyOptions::Read(const string& key, const string& msg)
{
    if (args.count(key) == 0)
        throw Invalid_option(msg);
                
    if (args[key] == "")
        throw Invalid_option("No value given for mandatory option: --" + key + "=???");
        
    // okay, option is valid.  Now remove it.
    string ret = args[key];
    args.erase(key);
    return ret;
}
        
bool EasyOptions::ReadBool(const string& key)
{
    if (args.count(key) == 0)
        return false;
        
    if (args[key] == "")
    {
        args.erase(key);
        return true;
    }
        
    throw Invalid_option(
        "Value should not be given for boolean option --" + key);
}
        
string EasyOptions::ReadWithDefault(const string& key, 
                                    const string& def)
{
    if (args.count(key) == 0)
        return def;
    if (args[key] == "")
        throw Invalid_option("Option requires a value: --" + key + 
            " (or omit, equivalent to --" + key + "=" + def);
    string ret = args[key];
    args.erase(key);
    return ret;
}


void EasyOptions::CheckEmpty()
{
    args.erase(""); // not worth complaining about this
    
    if (args.empty())
        return;
    
    string msg = "\nUnused arguments:\n" + stringify(*this);
    throw Invalid_option(msg);
}

ostream& operator<<(ostream& out, const EasyOptions& opts)
{
  for (map<string,string>::const_iterator i = opts.args.begin(); i != opts.args.end(); i++)
    {
      if (i->second == "")
        out << "--" << i->first << endl;
      else
        out << "--" << i->first << "='" << i->second << "'" << endl;
    }
  return out;
}

