# ddaParallelSubmission

****************************************************************************************************************************************
****************************************************************************************************************************************
Code for submitting parallel Discrete Dipole Approximation calculations to a remote computer.
****************************************************************************************************************************************
****************************************************************************************************************************************


A short description of each file:


****************************************************************************************************************************************
Directory -- root
****************************************************************************************************************************************

ddscatInputGen.py:

A python script adopted from a simliar script written by Zhongwei Hu designed to generate a parameter file for DDA from a simpler initialization.txt file and a shape.dat file. Can handle several parameter permutations--further details can be found in the ddscatInputGen.py file.



specDirGen.py:

A python script that generates all of the working directories (labeled WN, with N an integer 1 or larger) for all of the parallel DDA jobs to be run. Requires knowledge of the shape file, executable file, and scheduler submission script to be used in the calculation.



specRun.py:

A python script that runs all of the parllel jobs that are initialized by specDirGen.py.



specCollect.py:

A python script that collects the data from the parallel jobs run by specRun.py and aggregates it into one data file.



pIn.txt:

An example parameter initialization file to be used by ddscatInputGen.py.



