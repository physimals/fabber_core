/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

 Michael Chappell, FMRIB Image Analysis & IBME QuBIc groups

 Copyright (C) 2015 University of Oxford  */

/*  CCOPYRIGHT */

#ifndef __FABBER_CORE_H
#define __FABBER_CORE_H 1

#include "easylog.h"
#include "fabber_core.h"
#include "setup.h"
#include "utils.h"

#include <exception>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

/**
 * Run FABBER.
 *
 * @param argc Number of arguments.
 * @param argv Arguments controlling execution.
 * @return 0 if all went well, 1 otherwise.
 */
int execute(int argc, char** argv);

/**
 * Concatenate a vector of strings into a single string.
 * @param str_vector Vector of strings.
 * @param separator Separator for each string in the new string.
 * @return single string.
 */
string vectorToString(vector<string> str_vector, const char* separator = " ");

/**
 * Print usage information.
 */
void Usage();

/**
 * Print usage information for a specific forward model.
 * @param name Forward model name. If not known then an error will
 * be printed.
 */
void PrintFwdModelUsage(const string& name);

/**
 * Create a forward model using \ref FwdModelFactory.
 * @param name Forward model name.
 * @return pointer to forward model
 * @throws Invalid_option if the name is not known.
 */
FwdModel* GetFwdModel(const string& name);

/**
 * Create an inference technique using \ref InferenceTechniqueFactory.
 * @param name Inference technique name.
 * @return pointer to inference technqiue.
 * @throws Invalid_option if the name is not known.
 */
InferenceTechnique* GetInferenceTechnique(const string& name);

#endif // __FABBER_CORE_H
