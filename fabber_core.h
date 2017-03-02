/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

 Michael Chappell, FMRIB Image Analysis & IBME QuBIc groups

 Copyright (C) 2015 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

/**
 * Command line tool entry point for running Fabber.
 *
 * @param argc Number of arguments.
 * @param argv Command line arguments controlling execution.
 *
 * @return 0 if all went well, !=0 otherwise.
 */
int execute(int argc, char **argv);
