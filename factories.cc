/*  factories.cc - Assorted template factory, singleton factory and dispatcher class declarations.

     Mike Jackson, The University of Edinburgh & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

    Copyright (C) 2015 University of Oxford  */

#include "factories.h"

#include <vector>
#include <string>
#include <map>

using namespace std;

Dispatcher::Dispatcher()
{
    // No-op.
}

Dispatcher::~Dispatcher()
{
    functionMap_.clear();
}

void Dispatcher::Add(const string& name, Function function)
{
    functionMap_[name] = function;
}

void Dispatcher::Dispatch(const string& name)
{
    if (functionMap_.count(name)) {
        functionMap_[name]();
    }
}

vector<string> Dispatcher::GetNames()
{
    vector<string> names;
    for (map<string, Function>::iterator it = functionMap_.begin();
         it != functionMap_.end(); ++it) {
        names.push_back(it->first);
    }
    return names;
}

bool Dispatcher::HasName(const string& name)
{
    return (functionMap_.find(name) != functionMap_.end());
}
