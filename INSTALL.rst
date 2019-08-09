******************************
Compiling and installing DFTB+
******************************


Requirements
============

In order to compile DFTB+, you need the following software components:

* A Fortran 2003 compliant compiler

* A C-compiler

* CMake (version 3.5 or newer)

* GNU make

* LAPACK/BLAS libraries (or compatible equivalents)

Additionally there are optional requirements for some DFTB+ features:

* ScaLAPACK (version 2.0 or later) and an MPI aware Fortran compiler, if you
  want to build the MPI-parallelised version of the code

* The M4 preprocessor, if you want to build the MPI-parallelised version of the
  code

* The ARPACK or the ARPACK-ng library for excited state DFTB functionality

* The DftD3 dispersion library (if you need this dispersion model).

In order to execute the code tests and validate them against precalculated
results, you will additionally need:

* Python (version >= 2.6) with NumPy

* The Slater-Koster data used in the tests (see below)


Obtaining the source
====================

The source code of the last stable release can be downloaded from the `DFTB+
homepage <http://www.dftbplus.org>`_.

Alternatively you can clone the `public git repository
<https://github.com/dftbplus/dftbplus>`_. The tagged revisions correspond to
stable releases, while the master branch contains the latest development
version. As the project uses git-submodules, those must be additionally
downloaded ::

  git clone https://github.com/dftbplus/dftbplus.git
  cd dftbplus
  git submodule update --init --recursive


Optional extra components
~~~~~~~~~~~~~~~~~~~~~~~~~

Some optional software components are not distributed with the DFTB+ source
code. If these are required, but are not already installed on your system, then
we recommend you download these components by using the `get_opt_externals`
utility, e.g.::

  ./utils/get_opt_externals

This will download all license compatible optional external components. These
include the Slater-Koster (slako) data for testing the compiled code.

If you also wish to download and use any of the optional components which have
*conflicting licenses* (e.g. the `DftD3 library
<https://github.com/aradi/dftd3-lib>`_), you must explicitly request it::

  ./utils/get_opt_externals ALL

This will then prompt for confirmation when downloading components with other
licenses.

*Note*: if you include components with conflicting licenses into your
compilation of DFTB+, you are only allowed to use the resulting binary for your
personal research and are not permitted to distribute it.

For more information see the detailed help for this tool by issuing
``./utils/get_opt_externals -h``.


Compiling
=========

* Look at the `config.cmake` file and customise the global build parameters if
  necessary. (If you are unsure, leave the defaults as they are.)

* Customise the ``ARCH`` variable at the bottom of then `config.cmake` file. It
  should correspond to the prefix of one of the cmake config files in the `sys/`
  folder. Take the one which is closest to your system (and customise it
  according to your needs). Alternatively, you can also create your own file in
  the `sys/` folder.

  The build system will include the architecture dependent settings from the
  file `sys/${ARCH}.cmake` where ``${ARCH}`` is the value of the ``ARCH``
  variable you set in `config.cmake`.

* Create a build folder **outside** of the DFTB+ source tree and change to that
  that folder, e.g.::

    mkdir /tmp/build_dftb+
    cd /tmp/build_dftb+

* From the build folder Invoke CMake to configure the build. Pass the folder
  with the DFTB+ source code as argument, e.g. assuming the DFTB+ source is
  located in `~/dftbplus`, issue::

    cmake ~/dftbplus

* If the configuration was successful, invoke (from within the build folder)
  `make` to compile the code::

    make -j

  This will compile the code using several threads and showing only the most
  relevant information.

  If for debugging purposes you wish to see the exact compiling commands, you
  can execute a serial build with verbosity turned on, instead::

    make VERBOSE=1
  
* Note: The code can be compiled with distributed memory parallelism (MPI), but
  for smaller shared memory machines, you may find that the performance is
  better when using OpenMP parallelism only and an optimised thread aware BLAS
  library.


Testing DFTB+
=============

* After successful compilation, execute the code tests with ::

    make test

  You can also run the tests in parallel in order to speed this up.  If you use
  parallel testing, ensure that the number of OpenMP threads is reduced
  accordingly. As an example, assuming your workstation has 4 cores and you have
  set up the ``TEST_OMP_THREADS`` variable to ``2`` (in `config.cmake`), issue
  ::

    make test ARGS="-j2"

  for an OpenMP compiled binary running two tests simultaneously, each using 2
  cores.

  If you want to test the MPI enabled binary with more than one MPI-process, you
  should set the ``TEST_MPI_PROCS`` variable accordingly.

  Testing with hybrid (MPI/OpenMP) parallelism can be specified by setting both,
  the ``TEST_MPI_PROCS`` and ``TEST_OMP_THREADS`` variables, e.g::

    set(TEST_MPI_PROCS "2" CACHE STRING "Nr. of processes used for testing")
    set(TEST_OMP_THREADS "2" CACHE STRING "Nr. of OpeMP-threads used for testing")

  Note that efficient production use of the code in this mode may require
  process affinity (settings will depend on your specific MPI implementation).

  The ``TEST_MPI_PROCS`` and ``TEST_OMP_THREADS`` cache variables can also be
  updated dynamically just before starting the testing by invoking CMake with
  the appropriate ``-D`` options, e.g.::

    cmake -DTEST_MPI_PROCS=2 -DTEST_OMP_THREADS=2 ~/dftbplus
    make test


Installing DFTB+
================

* The compiled executables, libraries, module files etc. can be copied into an
  installation directory by ::

    make install

  where the destination directory can be configured by the variable
  ``CMAKE_INSTALL_PREFIX`` (in the `config.cmake` file). The default location is
  the `_install` subdirectory within the build directory.
