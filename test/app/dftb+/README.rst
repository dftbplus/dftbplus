********
AUTOTEST
********

The autotest suite contains some tests for the DFTB+ code with
precalculated results for automated comparison. In order to test your
DFTB+ binary on your system, you should execute the following command ::

  ./bin/autotest2 \
     -w <tempdir> \
     -p <binary> \
     -f <testfile> \
     -s A \
     -v

where:

* <tempdir> is an existing directory which should contain the results of the
            calculations to be validated (this directory is not deleted
            afterwards)

* <binary> is the path to your binary

* <testfile> contains the list of tests which should be calculated (you can use
             the 'tests' file in the this directory or create your own file)


If you have the DFTB+ source, you can let the build system execute the
tests for you ('ctest').


SLATER-KOSTER files
===================

You will have to provide Slater-Koster files for the autotest system. The cmake
system can download these directly ::

  mio-1-1
  hyb-0-2
  pbc-0-3
  rare-0-2



Important note
==============

The input files in the autotest system serve only to check if your binary works
correctly. The calculations made using them do not necessarily represent
meaningfull DFTB calculations for real physical systems. The geometries, the
tolerances, the spin constants etc. were often created in an arbitrary way
without any physical relevance. In some tests even the wrong Slater-Koster files
are used (e.g. scc SK-files for non-scc calculations). Therefore, DO NOT USE THE
INPUTS HERE AS TEMPLATES FOR YOUR INPUTS! If you need help in creating inputs,
consult the manual and the howtos on the dftb-plus.info webpage.


Adding additional tests
=======================

If you need to add additional autotests, perhaps as a result of
extending the functionality of DFTB+, there is a supplied script

bin/recalc_tests

which can generate the required _autotest.tag file. This is
traditionally generated using a DFTB+ binary compiled with the intel
compiler, and must be using a binary compiled for production use.

To use:

1. create an appropriate sub-directory inside `test/app/dftb+/` containing a
   `dftb_in.hsd` file for your test. This `dftb_in.hsd` file should be evaluated
   carefully for correctness, and MUST contain relevant parser version
   information.

2. The directory should also contain symbolic links to required Slater-Koster
   files from the appropriate external/slakos/origin/ directory.

3. in the directory test/app/dftb+ issue the command `./bin/recalc_tests
   ../../../_build/app/dftb+/dftb+ @tmp my_test` where my_test is the relative
   path from that directory to the new test directory.

4. Add an entry to the tests file in that directory.


For more complex test workflows, an optional `testrun.sh` file can be included
in the directory. This will be executed if present and can, for example, drive
multi-stage calculations or start required additional programs. In these cases,
the `recalc_tests` script will not be able to process these cases.
