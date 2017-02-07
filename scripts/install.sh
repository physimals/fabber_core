#!/bin/sh
#
# Build Fabber and install it
#
# Usage: install.sh [Debug|Release] [prefix]
#
# if prefix is not specified, uses $FSLDIR if set, $HOME if not

if [ -z $1 ] 
then
  TYPE=Debug
else
  TYPE=$1
fi

if [ -z $2 ] 
then
  if [ -z $FSLDIR ]
  then
    PREFIX=$HOME
  else
    PREFIX=$FSLDIR
  fi
else
  PREFIX=$2
fi

scriptdir=`dirname $0`
$scriptdir/build.sh $TYPE
cd build_$TYPE
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=$TYPE
make install


