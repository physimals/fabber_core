#!/bin/sh
#
# Build and install
#
# Usage: install.sh [debug|release] [prefix]
#
# if prefix is not specified, uses, in order of preference, $FABBERDIR, $FSLDIR
# or $HOME

if [ -z $1 ]
then
  TYPE=Debug
else
  TYPE=$1
fi

if [ -z $2 ]
then
  if [ -z $FABBERDIR ]
  then
    if [ -z $FSLDIR ]
    then
      PREFIX=$HOME
    else
      PREFIX=$FSLDIR
    fi
  else
    PREFIX=$FABBERDIR
  fi
else
  PREFIX=$2
fi

ORIGDIR=$PWD

scriptdir=`dirname $0`
$scriptdir/build.sh $TYPE
cd $scriptdir/../build_$TYPE
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=$TYPE
make install

cd $ORIGDIR
