#
# Global architecture independent build settings
#

#set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type (Release|RelWithDebInfo|Debug|MinSizeRel)")
# CMAKE_BUILD_TYPE is commented out in order to allow for multi-configuration builds. It will
# automatically default to RelWithDebInfo if used in a single configuration build. Uncomment or
# override it only if you want a non-default single configuration build.

option(WITH_OMP "Whether OpenMP thread parallisation should be enabled" TRUE)

option(WITH_MPI "Whether DFTB+ should support MPI-parallelism" FALSE)
# If you build an MPI-parallised binary, consider to set WITH_OMP (OpenMP thread parallelisaton) to
# FALSE unless you want hybrid parallelisation (for experts only).

option(WITH_ELSI "Whether DFTB+ with MPI-parallelism should use the ELSI libraries" FALSE)
# Works only with MPI-parallel build.

option(WITH_GPU "Whether DFTB+ should support GPU-acceleration via the MAGMA-library" FALSE)

option(WITH_TRANSPORT "Whether transport via libNEGF should be included." FALSE)
# Works only when building static libraries (see option BUILD_SHARED_LIBS)

option(WITH_SOCKETS "Whether socket communication should be allowed for" FALSE)

option(WITH_ARPACK "Whether the ARPACK library should be included (needed for TD-DFTB)" FALSE)
# Works only with non-MPI (serial) build, needed for Casida linear response

option(WITH_DFTD3 "Whether the DFTD3 library should be included" FALSE)
# NOTE: Due to the license of the DFTD3 library, the combined code must be distributed under the
# GPLv3 license (as opposed to the LGPLv3 license of the DFTB+ package)

option(WITH_MBD "Whether DFTB+ should be built with many-body-dispersion support" FALSE)

option(WITH_PLUMED "Whether metadynamics via the PLUMED2 library should be allowed for" FALSE)

option(WITH_API "Whether public API should be included and the DFTB+ library installed" TRUE)
# Turn this on, if you want to use the DFTB+ library to integrate DFTB+ into other software
# packages. (Otherwise only a stripped down version of the library without the public API is built.)
# This will also install necessary include and module files and further libraries needed to link the
# DFTB+ library.

option(WITH_PYTHON "Whether the Python components of DFTB+ should be tested and installed" TRUE)

option(BUILD_SHARED_LIBS "Whether the libraries built should be shared" FALSE)
# Turn this on, if the DFTB+ library (and other compiled libraries) should be shared libraries and
# dynamically linked to their applications. This results in smaller applications, but the libraries
# must be present at run-time (and the correct LD_LIBRARY_PATH environment variable must be set, so
# that they can be found by the operating system). If you want use the DFTB+ library from other
# software packages (see WITH_API option above), they may also require a shared library (e.g.
# calling DFTB+ functions from Python or Julia).


#
# Test environment settings
#
set(TEST_MPI_PROCS "1" CACHE STRING "Nr. of MPI processes used for testing")

set(TEST_OMP_THREADS "1" CACHE STRING "Nr. of OpenMP-threads used for testing")

# Command line used to launch the test code.
# The escaped variables (\${VARIABLE}) will be substituted by the corresponding CMake variables.
if(WITH_MPI)
  set(TEST_RUNNER_TEMPLATE "env OMP_NUM_THREADS=\${TEST_OMP_THREADS} mpiexec -n \${TEST_MPI_PROCS}"
    CACHE STRING "How to run the tests")
else()
  set(TEST_RUNNER_TEMPLATE "env OMP_NUM_THREADS=\${TEST_OMP_THREADS}" CACHE STRING
    "How to run the tests")
  set(MODES_RUNNER_TEMPLATE "env OMP_NUM_THREADS=\${TEST_OMP_THREADS}" CACHE STRING
    "How to run the modes code for tests")
endif()


#
# Installation options
#
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_install" CACHE STRING
  "Directory to install the compiled code into")

set(INSTALL_INCLUDEDIR "dftbplus" CACHE PATH
  "Name of the project specific sub-folder within the install folder for include files")

set(INSTALL_MODULEDIR "${INSTALL_INCLUDEDIR}/modfiles" CACHE PATH
  "Installation directory for Fortran module files (within the install folder for include files)")

set(PKGCONFIG_LANGUAGE "Fortran" CACHE STRING
  "Compiler and Linker language to assume when creating the pkg-config export file (C or Fortran)")
# The pkg-config export file (lib/pkgconfig/dftbplus.pc) contains the compiler and linker options
# needed to link the DFTB+ library to an application. (It can be queried with the pkg-config tool.)
# Depending on the language setting ("C" or "Fortran") you would get the flags for the case of using
# that compiler for the linking.


#
# Advanced options (e.g. for developers and packagers)
#

#set(TOOLCHAIN "gnu" CACHE STRING "Prefix of the toolchain file to be read from the sys/ folder")
# Uncomment and set it if you want to override the automatic, compiler based toolchain file
# selection.


set(HYBRID_CONFIG_METHODS "Submodule;Find;Fetch" CACHE STRING
  "Configuration methods to try in order to satisfy hybrid dependencies")
#
# This list can be used to control how hybrid dependencies (external dependencies which can
# optionally be built during the build process) are configured. The listed methods are applied in
# the following order:
#
# Submodule: Use the source in external/*/origin and build the dependency as part of the build
#     process. If the source is not present, try to retrieve it via the 'git submodule' command
#     (provided the source tree is a git repository and git is available)
#
# Find: Find the dependency as an already installed package in the system.
#
# Fetch: Fetch the source into the build folder and build the dependency as part of the build
#     process (works also in cases where the source tree is not a Git repository)
