#!/bin/bash
###   Compile everything
###   gfortran should be used
###
###   compile xdrfile for trr / xtc trajectories support
###
cd lib/xdrfile-1.1.4
./configure 
make clean
make 
cp src/.libs/libxdrfile.a ..
### Compiling analysis routines
cd ../../src 
make
echo "Done"

