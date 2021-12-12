*********************************************
Test framework for the modes vibrational code
*********************************************

AUTOTEST
========

The autotest suite contains some tests for the modes code from the DFTB+
project, with precalculated results for automated comparison. In order to test
this binary on your system, you should execute the following command:

.. code:: bash

  ./bin/autotest2 \
     -w <tempdir> \
     -p <binary> \
     -f <testfile> \
     -s A \
     -v

where

tempdir
             is an existing directory which should contain the results of the
             calculations to be validated (this directory is not deleted
             afterwards)

binary
             is the path to your binary

testfile
             contains the list of tests which should be calculated (you can use
             the 'tests' file in the this directory or create your own file)

If you have the DFTB+ source, you can let the cmake generated make system
execute the tests for you (`ctest -R modes` in the build directory).


SLATER-KOSTER files:
====================

You will have to provide Slater-Koster files for the autotest system (modes
requires masses for the atoms, which are currently read from these files unless
over-ridden). Please download the following SK-sets and extract them in the
external/slakos/origin/ subdirectory:

  mio-1-1

Note that if the `get_opt_externals` script has been run, it can find and
download this to the correct location.


Important note:
===============

The input files in the autotest system serve only to check if your binary works
correctly. The calculations made using them do not necessarily represent
meaningfull calculations for real physical systems.
