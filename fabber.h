
/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007-2012 University of Oxford  */

/*  CCOPYRIGHT */

// This header file is only needed if you want to link to fabber as a library.

#include <ostream>
#include <string>
#include <map>
#include <newmat.h>

using namespace std;
using namespace NEWMAT;

void fabber_library(const map<string,string>& argsIn,
                    const map<string,const Matrix*>& dataIn,
                    map<string,Matrix>& dataOut,
                    ostream& logOut);

