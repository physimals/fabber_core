#!/bin/sh
#
# Build distribution package for Fabber
#
# Usage: package.sh
#
# Prefix defaults to $FSLDIR if set, or $HOME if not

scriptdir=`dirname $0`
$scriptdir/build.sh release
cd $scriptdir/../build_release
cmake .. -DCMAKE_BUILD_TYPE=release
make package

scriptdir=`dirname $0`
$scriptdir/build.sh debug
cd $scriptdir/../build_debug
cmake .. -DCMAKE_BUILD_TYPE=debug
make package
