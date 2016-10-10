#!/bin/sh
#
# Build using Cmake

cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME
make install

