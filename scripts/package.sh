#!/bin/sh
#
# Build distribution package for Fabber
#
# Usage: package.sh 
#
# Prefix defaults to $FSLDIR if set, or $HOME if not

scriptdir=`dirname $0`
$scriptdir/build.sh Release
cd $scriptdir/../build_Release
cmake .. -DCMAKE_BUILD_TYPE=Release
make package


