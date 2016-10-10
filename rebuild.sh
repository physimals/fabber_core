#!/bin/sh
#
# Build using Cmake

rm -rf build
mkdir build
cd build
cmake ..
make all

