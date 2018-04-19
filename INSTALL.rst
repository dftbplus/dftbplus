******************************
Compiling and installing DFTB+
******************************


Requirements
============

In order to compile DFTB+, you need the following software components:

* A Fortran 2003 compliant compiler

* A C-compiler

* GNU make (version >= 3.79.1)

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

* Look at the makefiles in the `sys/` folder and find the one closest to your
  system. The suffix of the makefiles indicate the architecture, operating
  system and the compiler they have been written for. Copy the most suitable
  makefile to ``make.arch`` in the root directory of the source tree, e.g.::

      cp sys/make.x86_64-linux-gnu make.arch

* Adjust the settings in `make.arch` according to your system. Note that there
  are often separate settings in `make.arch` for compiling with and without
  MPI. The code is also usually compiled with openMP enabled.

* Open the file `make.config` and check the configuration options set there. In
  this file binary choices are defined as either 0 (false) or 1 (true).

* Build the binaries by issuing ::

     make

  in the root directory of your source tree. DFTB+ can be built in parallel, so
  you may use the ``-j`` option of `make` to specify the number of parallel
  build processes, e.g.::

    make -j4

  The build takes place in a separate directory `_build`. You can customise the
  build directory location within the `make.config` file (variable
  ``BUILDDIR``), but the default location is inside the root of the source tree.

* The code can be compiled with distributed memory parallelism (MPI), but for
  smaller shared memory machines, you may find that the performance is better
  when using OpenMP parallelism only and an optimised thread aware BLAS library.


Testing DFTB+
=============

* After successful compilation, execute the code tests with ::

    make test

  You can also run the tests in parallel (option ``-j``) in order to speed this
  up.  If you use parallel testing, ensure that the number of OpenMP threads is
  reduced accordingly. As an example, assuming your workstation has 4 cores, you
  could use::

    make -j2 test TEST_OMP_THREADS=2

  for an OpenMP compiled binary running two tests simultaneously, each using 2
  cores.

  If you want to test the MPI enabled binary with more than one MPI-process, you
  can set the TEST_MPI_PROCS variable accordingly e.g::

    make test TEST_MPI_PROCS=2

  Testing with hybrid (MPI/OpenMP) parallelism can be specified by setting both,
  the ``TEST_MPI_PROCS`` and ``TEST_OMP_THREADS`` variables, e.g::

    make test TEST_MPI_PROCS=2 TEST_OMP_THREADS=2

  Note that efficient production use of the code in this mode may require
  process affinity (settings will depend on your specific MPI implementation).

* The compiled executables can be copied into an installation directory by ::

    make install

  where the destination directory can be configured in the `make.config` file
  (set by the variable ``INSTALLDIR``).
