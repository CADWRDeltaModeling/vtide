VTide: Versatile Analysis of Tides
==================================

This code is a reworking of the Fisheries and Oceans Canada code described in Versatile Harmonic Tidal Analysis: Improvements and Applications by Foreman, Charniawsky and Ballantine. The goals for this project were:
* Facilitate prediction: the original versatile code required recourse to an older package to do prediction
* Open source solver: the new code links to LAPACK
* API amenable to python bindings: The api part is done, including separation of I/O and computation. Python bindings in progress. 

The new API may seem a bit non-F90, in the sense that the analysis routine does not pass any data via the backdoor route using module variables and without consolidating the analysis request and outputs into F90 derived types. The reason for the former is to allow multiple calls from different python instances. The reason for not using derived types is that they are not well supported by f2py.



Example command line
--------------------
cmake ..\src -DCMAKE_Fortran_COMPILER=ifort -G"Visual Studio 9 2008 Win64"


