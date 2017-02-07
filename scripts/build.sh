#!/bin/sh
#
# Make a clean build of Fabber
#
# Usage: build.sh [Debug|Release]

if [ -z $1 ] 
then
  TYPE=Debug
else
  TYPE=$1
fi

scriptdir=`dirname $0`
rm -rf $scriptdir/../build_$TYPE
mkdir $scriptdir/../build_$TYPE
cd $scriptdir/../build_$TYPE
cmake .. -DCMAKE_BUILD_TYPE=$TYPE
make

