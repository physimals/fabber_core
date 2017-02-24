/*  utils.cc - Assorted template factory, singleton factory and dispatcher class declarations.

     Mike Jackson, The University of Edinburgh & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

    Copyright (C) 2015 University of Oxford  */


#include "utils.h"

#include <miscmaths/miscmaths.h>

using namespace std;

using NEWMAT::Matrix;

using MISCMATHS::read_vest;
using MISCMATHS::read_ascii_matrix;

namespace fabber {

   Matrix read_matrix_file(std::string filename)
   {
	   // Detect if file contains a VEST matrix or not
	   try {
		   return read_vest(filename);
	   }
	   catch(...) {
		   // Do not care why this failed, if it was 'file not found'
		   // we will discover that now we try to read it as ASCII
		   return read_ascii_matrix(filename);
	   }
   }
}

Dispatcher::Dispatcher() {
  // No-op.
}

Dispatcher::~Dispatcher() { 
  functionMap_.clear(); 
}

void Dispatcher::Add(const string& name, Function function) {
  functionMap_[name] = function;
}

void Dispatcher::Dispatch(const string& name) {
  if (functionMap_.count(name)) {
    functionMap_[name]();
  }
}

vector<string> Dispatcher::GetNames() {
  vector<string> names;
  for (map<string, Function>::iterator it = functionMap_.begin(); 
       it != functionMap_.end(); it++) {
    names.push_back(it->first);
  }
  return names;
}

bool Dispatcher::HasName(const string& name) {
  return (functionMap_.find(name) != functionMap_.end());
}
