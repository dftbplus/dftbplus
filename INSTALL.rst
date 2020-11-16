*****************************
Building and installing DFTB+
*****************************

If you have problems with the build, you can find suggestions for some
frequently occuring scenarios in the `Troubleshooting <#troubleshooting>`_
section at the bottom.


Requirements
============

In order to compile DFTB+, you need the following software components:

* A Fortran 2003 compliant compiler

* A C-compiler

* CMake (version 3.16 or newer)

* GNU make

* LAPACK/BLAS libraries (or compatible equivalents)

* Python (version >= 3.2) for the source preprocessor


Optional extra dependencies
---------------------------

Additionally there are optional requirements for some DFTB+ features:

* ScaLAPACK (version 2.0 or later) and a Fortran aware MPI framework, if you
  want to build the MPI-parallelised version of the code.

* In addition to ScaLAPACK, it is recommended to use the `ELSI
  <https://wordpress.elsi-interchange.org/>`_ library for large scale systems
  (version 2.6.x of the library, with partial support of 2.5.0). If ELSI was
  compiled with PEXSI included, you will also need a C++ compiler.

* The ARPACK-ng library if using the excited state DFTB functionality.

* The `MAGMA <http://icl.cs.utk.edu/magma/>`_ library for GPU accelerated
  computation.

* The `PLUMED2 <https://github.com/plumed/plumed2>`_ library for metadynamics
  simulations. If you build DFTB+ with MPI, the linked PLUMED library must be
  also MPI-aware (and must have been built with the same MPI-framework as
  DFTB+).


External library requirements
-----------------------------

* **Make sure that all external libraries are compiled with the same kind models
  for the numeric variables** (same integer size and floating point precision)
  as DFTB+. Also, they should preferably have been built with the same compiler
  and with similar compiler flags to DFTB+. (See the Troubleshooting section for
  further information.)

* External libraries in non-standard locations (as is typical on many
  HPC-systems using environment modules) can only be reliable found by CMake if
  their library path occurs in the ``CMAKE_PREFIX_PATH`` environment
  variable. **Make sure that your CMAKE_PREFIX_PATH environment variable
  contains all relevant library paths!**


Requirements for testing DFTB+
------------------------------

In order to execute the code tests and validate them against precalculated
results, you will additionally need:

* Python (version >= 3.2) with NumPy

* The Slater-Koster data used in the tests (see below)


Tested build environments
-------------------------

DFTB+ is regularly built and tested for both serial and MPI environments on the
following architectures:

+---------------+----------------------+-------------+------------------+-----+
| Architecture  | Compiler             | MPI         | Ext. libraries   |Notes|
+===============+======================+=============+==================+=====+
| x86_64 /      | GNU Fortran/C 7.5    | OpenMPI 2.1 | OpenBlas 0.3.7,  |     |
| Linux         |                      |             | ScaLAPACK 2.1    |     |
|               |                      |             | ELSI 2.6.1       |     |
+---------------+----------------------+-------------+------------------+-----+
| x86_64 /      | GNU Fortran/C 10.1   | OpenMPI 4.0 | OpenBlas 0.3.10, |     |
| Linux         |                      |             | ScaLAPACK 2.1    |     |
|               |                      |             | ELSI 2.6.1       |     |
+---------------+----------------------+-------------+------------------+-----+
| x86_64 /      | Intel Fortran/C 18.0 | MPICH 3.2   | MKL 18.0         |     |
| Linux         |                      |             | ELSI 2.6.1       |     |
+---------------+----------------------+-------------+------------------+-----+
| x86_64 /      | Intel Fortran/C 19.0 | MPICH 3.3   | MKL 19.0         |     |
| Linux         |                      |             | ELSI 2.6.1       |     |
+---------------+----------------------+-------------+------------------+-----+
| x86_64 /      | NAG Fortran 7.0      | MPICH 3.3   | OpenBlas 0.3.7   |     |
| Linux         | GNU C 9.2            |             | ScaLAPACK 2.1    |     |
|               |                      |             | ELSI 2.5.0       |     |
+---------------+----------------------+-------------+------------------+-----+
| x86_64 /      | GNU Fortran/C 8.4    | --          | OpenBlas 0.3.10  | [1] |
| OS X          |                      |             |                  |     |
+---------------+----------------------+-------------+------------------+-----+

All builds are also tested with the optional ARPACK-NG 3.7 and PLUMED 2.5
libraries.

Notes:

[1] Only partial testing of the serial version.


Obtaining the source
====================

The source code of the last stable release can be downloaded from the `DFTB+
homepage <https://www.dftbplus.org/download/dftb-stable/>`_.

Alternatively you can clone the `public git repository
<https://github.com/dftbplus/dftbplus>`_. The tagged revisions correspond to
stable releases, while the default branch contains the latest development
version. ::

  git clone https://github.com/dftbplus/dftbplus.git
  cd dftbplus

The project uses git-submodules for some external dependencies, which will be
automatically retrieved during configuration.


Optional extra components
-------------------------

Some optional software components are not distributed with the DFTB+ source code
and are also not retrieved automatically. If these are required, you can
download these components by using the `get_opt_externals` utility, e.g.::

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


Building
========

**Important note:** CMake caches its variables in the `CMakeCache.txt` file in
the build folder (e.g. ``_build/CMakeCache.txt``). Make sure to delete this file
before re-running CMake if you have changed variables in `config.cmake` or in
the toolchain files in the `sys/` folder. (Deleting the `CMakeCache.txt` file is
not necessary if you change a variable via the ``-D`` command line option.)

In order to build DFTB+ carry out the following steps:

* Inspect the `config.cmake` file and customise the global build parameters. (If
  you are unsure, leave the defaults as they are.)

* Invoke CMake to configure the build. Specify the installation destination
  (e.g. ``$HOME/opt/dftb+``) and pass an arbitrary folder (e.g. ``_build``) for
  the build and the directory containing the source files (e.g. ``.``) as
  arguments to CMake. Additionally define your Fortran and C compilers as
  environment variables, e.g. (in a BASH compatible shell)::

    FC=gfortran CC=gcc cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt/dftb+ -B _build .

  Based on the detected compilers, the build system will read further settings
  from a corresponding toolchain file in the `sys/` folder. Either from a
  compiler specific one (e.g. `gnu.cmake`, `intel.cmake`, etc.) or the generic
  one (`generic.cmake`) if the detected compiler combination does not correspond
  to any of the specific settings. The selected toolchain is indicated in the
  CMake output. (The toolchain file selection can be manually overriden by
  setting the ``TOOLCHAIN`` CMake variable.)

  You may adjust any CMake variable defined in `config.make` or in the
  toolchain files by either modifying the files directly or by setting
  (overriding) the variable via the ``-D`` command line option. For example, in
  order to use the MKL-library with the GNU-compiler, you would have to override
  the ``LAPACK_LIBRARY`` variable with the CMake command line argument ``-D``::

    -DLAPACK_LIBRARY="mkl_gf_lp64;mkl_gnu_thread;mkl_core"

  When needed, you can specify the complete path to a libray or pass linker
  options as defined variables, e.g.::

    -DLAPACK_LIBRARY="/opt/openblas/libopenblas.a"
    -DLAPACK_LIBRARY="-Wl,--start-group -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group"

  By default CMake searches for the external libraries in the paths specified in
  the ``CMAKE_PREFIX_PATH`` environment variable. **Make sure that your
  CMAKE_PREFIX_PATH environment variable is set up correctly and contains
  all the relevant paths** when configuring the project, e.g. ::

    CMAKE_PREFIX_PATH=/opt/elsi:/opt/custom-openblas cmake [...] -B _build .

  Some of the external library finders also offer special ``_LIBRARY_DIR`` CMake
  variables for setting search paths, e.g. ::

    -DLAPACK_LIBRARY_DIR=/opt/custom-openblas

  Setting those variables is not normally necessary, if the right search path is
  already present in the ``CMAKE_PREFIX_PATH`` environment variable.


* If the configuration was successful, start the build by ::

    cmake --build _build -- -j

  This will compile the code using several threads and showing only the most
  relevant information.

  If, for debugging purposes, you wish to see the exact compiling commands, you
  should execute a serial build with verbosity turned on instead::

    cmake --build _build -- VERBOSE=1

* Note: The code can be compiled with distributed memory parallelism (MPI), but
  for smaller shared memory machines, you may find that the performance is
  better when using OpenMP parallelism only and an optimised thread aware BLAS
  library is used.


Testing DFTB+
=============

* After successful compilation, change to the build folder and execute the code
  tests::

    pushd _build
    ctest
    popd

  You can also run the tests in parallel in order to speed this up.  If you use
  parallel testing, ensure that the number of OpenMP threads is reduced
  accordingly. As an example, assuming your workstation has 4 cores and you have
  set up the ``TEST_OMP_THREADS`` variable to ``2`` (in `config.cmake`), issue
  ::

    ctest -j2

  for an OpenMP compiled binary running two tests simultaneously, each using 2
  cores.

  If you want to test the MPI enabled binary with more than one MPI-process, you
  should set the ``TEST_MPI_PROCS`` variable accordingly.

  Testing with hybrid (MPI/OpenMP) parallelism can be specified by setting both,
  the ``TEST_MPI_PROCS`` and ``TEST_OMP_THREADS`` variables, e.g::

    set(TEST_MPI_PROCS "2" CACHE STRING "Nr. of processes used for testing")
    set(TEST_OMP_THREADS "2" CACHE STRING "Nr. of OMP-threads used for testing")

  Note that efficient production use of the code in this mode may require
  process affinity (settings will depend on your specific MPI implementation).

  The ``TEST_MPI_PROCS`` and ``TEST_OMP_THREADS`` cache variables can be updated
  or changed also after the compilation by invoking CMake with the appropriate
  ``-D`` options, e.g.::

    cmake -B _build -DTEST_MPI_PROCS=2 -DTEST_OMP_THREADS=2 .
    pushd _build; ctest; popd


Installing DFTB+
================

* The compiled executables, libraries, module files etc. can be copied into an
  installation directory by ::

    cmake --install _build

  where the destination directory can be configured by the variable
  ``CMAKE_INSTALL_PREFIX`` (in the `config.cmake` file). The default location is
  the `_install` subdirectory within the build directory.


Using DFTB+ as a library
========================

DFTB+ can be also be used as a library and linked into other simulation software
packages. In order to compile the library with its public API, make sure to set
the ``WITH_API`` option to ``TRUE`` in the CMake config file
`config.cmake`. When you install the program, it will also install the DFTB+
library, the C-include file and the Fortran module files, which are necessary
for linking DFTB+ with C and Fortran programs.


Linking the library in CMake based builds
-----------------------------------------

This is the prefered way of invoking the DFTB+ library into your project.  In
CMake based projects you can directly use the CMake export file of DFTB+, which
is installed in the `lib/cmake/dftbplus/` folder in the installation folder. It
exports the target ``DftbPlus::DftbPlus`` which you can use to obtain all
necessary compiler, include and linking options. Your projects `CMakeLists.txt`,
should like something like below::

  project(DftbPlusTest LANGUAGES Fortran C)
  find_package(DftbPlus REQUIRED)
  add_executable(testprogram testprogram.f90)
  target_link(testprogram DftbPlus::DftbPlus)

Note, that this will link all libraries in the correct order, which where
compiled during the DFTB+ build (e.g. libdftd3, libnegf, etc.). It will
additionally contain target dependencies on the external libraries needed to
create standalone applications with DFTB+ (e.g. ``LAPACK::LAPACK``,
``Scalapack::Scalapack``, ``Arpack::Arpack``, ``Plumed::Plumed``,
``Magma::Magma``, etc.). You can either use the CMake find-modules shipped with
the DFTB+ source to find those libraries (and to define the corresponding
targets) or create your own ones, provided they define the appropriate CMake
targets. The ELSI library offers a CMake export file providing the
``elsi::elsi`` target. Make sure, that CMake can find this export file if the
DFTB+ library was compiled with ELSI support (e.g. by setting up the environment
variable ``CMAKE_PREFIX_PATH`` correctly).


Linking the library in non-CMake based builds
---------------------------------------------

Depending on the choice of external components and whether you want to link
DFTB+ to a C or a Fortran binary, you may need different compilation flags and
linker options. You can look up the necessary compiler flags and linker options
in the `dftbplus.pc` pkg-config file, which is usually installed into the
`lib/pkgconfig` folder in the installation directory. You can either inspect the
file directly, or use the ``pkg-config`` tool::

  export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:DFTBPLUS_INSTALL_FOLDER/lib/pkgconfig
  pkg-config --cflags dftbplus   # compilation flags (e.g. include options)
  pkg-config --libs dftbplus     # library linking options
  pkg-config --static --libs dftbplus   # library linking options for static linking

Note, that the flags and libraries shown are either for linking with Fortran or
with C, depending on the value of the configuration option
``PKGCONFIG_LANGUAGE``.

If you compile DFTB+ with ELSI, PLUMED or MAGMA-support, make sure that
pkg-config can also find their respective pkconfig files, as those libraries are
declared as dependencies in the DFTB+ pkg-config file. For external dependencies
without pkg-config files (e.g. mbd, negf) the options for linking those
libraries can not be queried via pkg-config and must be added manually.


Generating developer documentation
==================================

Developer documentation can be generated using the FORD source code
documentation generator by issuing ::

  cd doc/dftb+/ford && ford dftbplus-project-file.md

in the main source directory. The documentation will be created in the
`doc/dftb+/ford/doc` folder.


Developer build instructions
============================

You should avoid customizing the build by directly changing variables in the
CMake config files, as your changes may accidently be checked in into the
repository. Instead, create a customized CMake config file, where you
pre-populate the appropriate cache variables. Then use the `-C` option to load
that file::

  FC=gfortran CC=gcc cmake -C custom.cmake -B _build .

The customized config file is read by CMake before the compiler detection
stage. If your config file contains toolchain dependent options, consider
defining the ``DFTBPPLUS_TOOLCHAIN`` environment variable and query it in your
config file.


Advanced build configuration (e.g. for packagers)
=================================================

Controlling the toolchain file selection
----------------------------------------

You can override the toolchain file, and select a different provided case,
passing the ``-DTOOLCHAIN`` option with the relevant name, e.g.::

  -DTOOLCHAIN=gnu

or by setting the toolchain name in the ``DFTBPLUS_TOOLCHAIN`` environment
variable. If you want to load an external toolchain file instead of one from the
source tree, you can specify the file path with the ``-DTOOLCHAIN_FILE`` option
::

  -DTOOLCHAIN_FILE=/some/path/myintel.cmake

or with the ``DFTBPLUS_TOOLCHAIN_FILE`` environment variable.

Similarly, you can also use an alternative build config file instead of
`config.cmake` in the source tree by specifying it with the
``-DBUILD_CONFIG_FILE`` option or by defining the ``DFTBPLUS_BUILD_CONFIG_FILE``
environment variable.


Preventing the download of external sources
-------------------------------------------

Depending on the value of the ``HYBRID_CONFIG_METHODS`` configuration variable,
some dependencies (e.g. mbd, negf, mpifx, scalapackfx) are automatically
downloaded during the configuration phase and built during the DFTB+ build
process. If you want to ensure that nothing gets downloaded during the build,
pass the variable definition ::

  -DHYBRID_CONFIG_METHODS="Find"

to CMake during the configuration. In this case, CMake will only try to find
those dependencies on the system (by searching in the standard system paths and
in the locations defined in the environment variable ``CMAKE_PREFIX_PATH``) and
stop if some components were not found.


Troubleshooting
===============

* **CMake finds the wrong compiler**

  CMake should be guided with the help of the environment variables ``FC``,
  ``CC`` (and eventually ``CXX``) to make sure it uses the right compilers,
  e.g. ::

    FC=gfortran CC=gcc cmake [...]


* **CMake fails to find a library / finds the wrong version of a library**

  In most cases this is due to a misconfigured ``CMAKE_PREFIX_PATH`` environment
  variable. It is essential, that ``CMAKE_PREFIX_PATH`` contains all paths
  (besides default system paths), which CMake should search when trying to find
  a library. Extend the library path if needed, e.g. ::

    CMAKE_PREFIX_PATH="/opt/somelib:${CMAKE_PREFIX_PATH}" cmake [...]


* **ScaLAPACK detection on Ubuntu 20.4 LTS fails**

  The OpenMPI version of ScaLAPACK on Ubuntu 20.4 LTS exports an incorrect CMake
  config file (as of October 2020), which refers to an non-existent
  library. Instead, set the library name with the ``SCALAPACK_LIBRARY`` variable
  explicitely, e.g. ::

    cmake -DSCALAPACK_LIBRARY=scalapack-openmpi [...]

  which should fix the problem.


* **My library settings in a "_LIBRARIES" variable are ignored**

  In order to be consistent with the naming scheme suggested by the CMake
  documentation, all library related cache variables have been changed to
  singular nouns, e.g. ::

    cmake -DSCALAPACK_LIBRARY=scalapack-openmpi [...]

  **instead** of the previous ::

    cmake -DSCALAPACK_LIBRARIES=scalapack-openmpi [...]


* **Fortran libraries compiled with the Intel compiler can not be linked**

  In order to enforce compliance with the Fortran 2003 standard (e.g. allowing
  the automatic allocation of arrays in expressions), DFTB+ passes the
  ``-standard-semantics`` option to the Intel compiler. All external modern
  Fortran dependencies (e.g. ELSI) must also be compiled by using the
  ``-standard-semantics`` or the ``-assume realloc_lhs`` option to ensure
  correct linking.
