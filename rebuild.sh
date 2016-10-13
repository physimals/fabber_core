#!/bin/sh
#
# Build using Cmake

if [ -z $1 ] 
then
  TYPE=Debug
else
  TYPE=$1
fi

rm -rf $TYPE
mkdir $TYPE
cd $TYPE
cmake .. -DCMAKE_BUILD_TYPE=$TYPE
make all

