/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

 Michael Chappell, FMRIB Image Analysis & IBME QuBIc groups

 Copyright (C) 2015 University of Oxford  */

/*  CCOPYRIGHT */

#ifndef __FABBER_CORE_H
#define __FABBER_CORE_H 1

/**
 * Run FABBER.
 *
 * @param argc Number of arguments.
 * @param argv Arguments controlling execution.
 * @return 0 if all went well, 1 otherwise.
 */
int execute(int argc, char** argv);

#endif // __FABBER_CORE_H
