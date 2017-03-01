#!/bin/sh
#
# Build distribution packages
#
# Usage: package.sh

ORIGDIR=$PWD
scriptdir=`dirname $0`

$scriptdir/build.sh release
cd $scriptdir/../build_release
make package
cd ..

$scriptdir/build.sh debug
cd $scriptdir/../build_debug
make package
cd ..

cd $ORIGDIR
