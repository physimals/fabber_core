#!/bin/sh
#
# Install 

if [ -z $1 ] 
then
  TYPE=Debug
else
  TYPE=$1
fi

if [ -z $PREFIX ]
then
  PREFIX=$HOME
fi

cd $TYPE
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX
make install

