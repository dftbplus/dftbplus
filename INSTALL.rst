******************************
Compiling and installing DFTB+
******************************


Requirements
============

In order to compile DFTB+, you need the following software components:

* A Fortran 2003 compliant compiler

* A C-compiler (required if building with the socket interface enabled or if C
  language API bindings are required)

* CMake (version 3.5 or newer)

* GNU make

* LAPACK/BLAS libraries (or compatible equivalents)

* Python (version >= 3.2) for the source preprocessor

Additionally there are optional requirements for some DFTB+ features:

* ScaLAPACK (version 2.0 or later) and an MPI aware Fortran compiler, if you
  want to build the MPI-parallelised version of the code

* In addition to ScaLAPACK, the `ELSI
  <https://wordpress.elsi-interchange.org/>`_ library for large scale systems
  can optionally also be used (version 2.5.0 or later of the library).

* The ARPACK or the ARPACK-ng library for excited state DFTB functionality

* The `MAGMA <http://icl.cs.utk.edu/magma/>`_ library for GPU accelerated
  computation.

* The `PLUMED2 <https://github.com/plumed/plumed2>` library for metadynamics
  simulations. If you build DFTB+ with MPI, the linked PLUMED library must be
  also MPI-aware (and must have been built with the same MPI-framework as
  DFTB+).

For external libraries, make sure that they are compiled with the same precision
models for the variables (same integer and floating point values).

In order to execute the code tests and validate them against precalculated
results, you will additionally need:

* Python (version >= 3.2) with NumPy

* The Slater-Koster data used in the tests (see below)


Obtaining the source
====================

The source code of the last stable release can be downloaded from the `DFTB+
homepage <http://www.dftbplus.org>`_.

Alternatively you can clone the `public git repository
<https://github.com/dftbplus/dftbplus>`_. The tagged revisions correspond to
stable releases, while the default branch contains the latest development
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


Building
========

* Look at the `config.cmake` file and customise the global build parameters if
  necessary. (If you are unsure, leave the defaults as they are.)

* Create a build folder (e.g. `_build`) either in the DFTB+ source tree or
  somewhere else outside of it and change to that folder, e.g.::

    mkdir _build
    cd _build

* From the build folder invoke CMake to configure the build. In order to ensure,
  that CMake chooses the right compilers, pass the Fortran and C compilers as
  environment variables, e.g.::

    env FC=gfortran CC=gcc cmake ..

  If you want to compile the code with MPI-support, you should pass the name of
  the mpi compiler-wrapper as Fortran compiler, e.g.::

    env FC=mpifort CC=gcc cmake ..

  Based on the detected compilers, the build system will read further settings
  from the corresponding pre-configured toolchain file in the `sys/` folder.

  You may adjust any variables in the `config.make` file and the toolchain file
  by either changing them directly in the files or by overriding them via the
  ``-D`` command line option. For example, in order to change the name of the
  LAPACK library, you could override the ``LAPACK_LIBRARIES`` variable with::

    env FC=gfortran CC=gcc cmake -DLAPACK_LIBRARIES=openblas-custom ..

  CMake automatically searches for the external libraries in the paths specified
  in the ``CMAKE_PREFIX_PATH`` environment variable. Make sure that it is set up
  correctly in your build environment. Alternatively, the DFTB+ build system
  offers a ``*_LIBRARY_DIRS`` variable for each external library, where you can
  add path hints for the library search. For example ::

    env FC=gfortran CC=gcc cmake -DLAPACK_LIBRARIES=openblas \
        -DLAPACK_LIBRARY_DIRS=/opt/openblas/lib

  would advise CMake to look also in the `/opt/openblas/lib` directory for the
  LAPACK (OpenBlas) library.

  If you want to select a different toolchain as choosen by the build system,
  you can override the ``TOOLCHAIN`` variable::

    env FC=ifort CC=gcc cmake -DTOOLCHAIN=generic ..

  Alternatively, you can also specify a custom toolchain file at an arbitraty
  location with the ``TOOLCHAIN_FILE`` variable, e.g.::

    env FC=ifort CC=gcc cmake -DTOOLCHAIN_FILE=/somepath/myintelgnu.cmake ..

  You can also override the ``TOOLCHAIN`` and ``TOOLCHAIN_FILE`` variables by
  defining the environment variables ``DFTBPLUS_TOOLCHAIN`` or
  ``DFTBPLUS_TOOCHAIN_FILE``.


* If the configuration was successful, invoke (from within the build folder)
  `make` to compile the code::

    make -j

  This will compile the code using several threads and showing only the most
  relevant information.

  If, for debugging purposes, you wish to see the exact compiling commands, you
  should execute a serial build with verbosity turned on instead::

    make VERBOSE=1
  
* Note: The code can be compiled with distributed memory parallelism (MPI), but
  for smaller shared memory machines, you may find that the performance is
  better when using OpenMP parallelism only and an optimised thread aware BLAS
  library.


Testing DFTB+
=============

* After successful compilation, execute the code tests with ::

    ctest

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

    cmake -DTEST_MPI_PROCS=2 -DTEST_OMP_THREADS=2 ..
    ctest


Installing DFTB+
================

* The compiled executables, libraries, module files etc. can be copied into an
  installation directory by ::

    make install

  where the destination directory can be configured by the variable
  ``CMAKE_INSTALL_PREFIX`` (in the `config.cmake` file). The default location is
  the `_install` subdirectory within the build directory.



Using DFTB+ as a library
========================

DFTB+ can be also used as a library and linked with other simulation software
packages. In order to compile the library with the public API, make sure to set
the ``WITH_API`` option to ``TRUE`` in the CMake config file
`config.cmake`. When you install the program, it will also install the DFTB+
library (`libdftbplus.a`), the C-include file and the Fortran module files,
which are necessary for linking DFTB+ with C and Fortran programs.


Linking the library in non-CMake based builds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

If you compile DFTB+ with ELSI-support, make sure that pkg-config can find
ELSI's own pkgconfig file, as it is declared as dependency in the DFTB+
pkg-config file.


Linking the library in CMake based builds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you use CMake to build your project, you can directly use the CMake
configuration file installed by DFTB+ into the `lib/cmake/DftbPlus/` folder in
the installation root directory. It exports the target ``DftbPlus::dftbplus``
which you can use to obtain compiler and linking options. For example, in your
projects `CMakeLists.txt`, you could have something like::

  project(dftbplus_libtest LANGUAGES Fortran C)
  find_package(DftbPlus REQUIRED)
  add_executable(testprogram testprogram.f90)
  target_link(testprogram DftbPlus::dftbplus)

Note, that this will link all libraries in the correct order, which where
compiled during the DFTB+ build (e.g. libdftd3, libnegf, etc.). It will also
contain the link dependencies on the external libraries needed to create
standalone applications with DFTB+ (e.g. lapack, scalapack). You must make sure,
that CMake can find those libraries, when linking the
application. Alternatively, you may use CMake to find them at the locations,
where they were found during the DFTB+ build. The variables
``DftbPlus_EXTERNAL_LIBRARIES`` and ``DftbPlus_EXTERNAL_LIBRARY_DIRS`` contain
all external libraries and the directories, where they have been found. In order
to make sure, CMake finds them, you could turn them into targets in your CMake::

  project(dftbplus_libtest LANGUAGES Fortran)

  find_package(DftbPlus REQUIRED)

  foreach(lib IN LISTS DftbPlus_EXTERNAL_LIBRARIES)
    find_library(LIBPATH ${lib} HINTS ${DftbPlus_EXTERNAL_LIBRARY_DIRS})
    if(LIBPATH)
      message(STATUS "Found library ${LIBPATH}")
      add_library(${lib} IMPORTED UNKNOWN)
      set_target_properties(${lib} PROPERTIES IMPORTED_LOCATION ${LIBPATH})
    else()
      message(FATAL_ERROR
        "Could not find library '${lib}' using library path hints '${libpaths}'")
    endif()
    unset(LIBPATH CACHE)
  endforeach()

  add_executable(testprogram testprogram.f90)
  target_link_libraries(testprogram DftbPlus::dftbplus)

If you compile DFTB+ with ELSI support, make sure that CMake can find ELSI's own
CMake configuration file, as it is declared as dependency in the DFTB+ Cmake
config file.


Generating developer documentation
==================================

Developer documentation can be generated using the FORD source code
documentation generator by issuing ::

  cd doc/dftb+/ford && ford dftbplus-project-file.md

in the main source directory. The documentation will be created in the
`doc/dftb+/ford/doc` folder.
