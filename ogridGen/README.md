This directory contains FORTRAN scripts that will generate an O-grid using the output from the Bezier smoother.  It assumes the input file is named "airfoil.p3d" and will generate a file "grid.in" which will be a three planed 2-D airfoil grid.  

The files "addendpoints.F90" and "addplanes.F90" need to be compiled with the following commands

`gfortran -fdefault-real-8 addendpoints.F90 -o addendpoints`

`gfortran -fdefault-real-8 addplanes.F90 -o addplanes`

The "hypgen.i" file should be edited for desired "k" direction properties.

The script "gengrid" assumes the executables "addendpoints" and "addplanes" are called locally but can be editted if they are in PATH
