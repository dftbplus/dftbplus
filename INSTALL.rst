******************************
Compiling and installing DFTB+
******************************


Requirements
============

In order to install DFTB+, you need the following software components:

* A Fortran 2003 compliant compiler

* A C-compiler

* GNU make (version >= 3.79.1)

* LAPACK/BLAS libraries (or compatible equivalents)

* Optionally: ScaLAPACK (version 2.0 or later), if you build the
  MPI-parallelised version of the code

* Optionally: the ARPACK or ARPACK-ng library for excited state DFTB
  functionality; the DftD3 dispersion library (if you need this dispersion
  model).

In order to execute the tests and compare them against precalculated results,
you will additionally need:

* Python (version >= 2.6) with NumPy


Obtaining the source
====================

The source code can be downloaded from the `DFTB+ homepage
<http://www.dftbplus.org>`_.

Alternatively you can clone the `public git repository
<https://github.com/dftbplus/dftbplus>`_. (The tagged revisions correspond to
stable releases, while the master branch contains the latest development
version.) As the project uses git-submodules, those must be additionally
downloaded ::

  git clone https://github.com/dftbplus/dftbplus.git
  cd dftbplus
  git submodule update --init --recursive

Some optional software components (e.g. the `DftD3 library
<https://github.com/aradi/dftd3-lib>`_) are not distributed with the DFTB+
source code, but can be included during the DFTB+ compilation if they are not
installed on your system. You can download those optional software components by
using the `get_opt_externals` utility, e.g.::

  ./utils/get_opt_externals dftd3

The Slater-Koster data needed for testing can also be downloaded by using
this tool::

  ./utils/get_opt_externals slakos

See detailed help for this tool by issuing ``./utils/get_opt_externals -h``.


Compile
=======

* Look at the makefiles in the `sys/` folder and find the one closest to your
  system. The suffix of the makefiles indicate the architecture, operating
  system and the compiler they have been written for. Copy the most suitable
  makefile to ``make.arch`` in the root directory of the source tree, e.g.::

      cp sys/make.x86_64-linux-gnu make.arch

* Adjust the settings in `make.arch` according to your system.

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

* After successful compilation, execute the tests with ::

    make test

  You may also run the tests in parallel (option ``-j``) in order to speed this
  up.  If you use parallel testing, ensure that the number of OpenMP threads is
  set to be 1 via the ``OMP_NUM_THREADS`` environment variable before starting
  the tests, e.g.::

    export OMP_NUM_THREADS=1

  if using the bash shell. If you want to test the MPI-binary with more than one
  processes, you can set the TESTPROC variable accordingly e.g::

    make test TESTPROC=2

* The compiled executables can be copied into an installation directory by ::

    make install

  where the destination directory can be configured in the `make.config` file
  (set by the variable ``INSTALLDIR``).
