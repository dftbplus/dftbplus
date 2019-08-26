#
# Global architecture independent build settings
#
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Release|Debug)")

set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_install" CACHE STRING
  "Directory to install the compiled code into")

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
# Works only with non-MPI (serial) build

option(WITH_DFTD3 "Whether the DFTD3 library should be included" FALSE)
# NOTE: Due to the license of the DFTD3 library, the combined code must be distributed under the
# GPLv3 license (as opposed to the LGPLv3 license of the DFTB+ package)

option(BUILD_SHARED_LIBS "Whether the libraries built should be shared" FALSE)
# Turn this on, if the DFTB+ library (and other compiled libraries) should be shared libraries and
# dynamically linked to their applications. This results in smaller applications, but the libraries
# must be present at run-time (and the correct LD_LIBRARY_PATH environment variable must be set, so
# that they can be found by the operating system). If you want use the DFTB+ library from other
# software packages (see BUILD_API option below), they may also require a shared library (e.g.
# calling DFTB+ functions from Python or Julia).

option(BUILD_API "Whether DFTB+ library with high-level API should be built and installed" FALSE)
# Turn this on, if you want to use the DFTB+ library to integrate DFTB+ into other software
# packages. (Otherwise only a stripped down version of the library without the public API is built.)
# This will also install necessary include and module files and further libraries needed to link the
# DFTB+ library.


#
# Architecture/compiler specific build settings
#
set(ARCH "x86_64-linux-gnu" CACHE STRING
  "Selects which architecture dependent settings should be used")

# Include compiler dependent build settings from the sys-directory
include(${CMAKE_SOURCE_DIR}/sys/${ARCH}.cmake)


# Test environment settings
set(TEST_MPI_PROCS "1" CACHE STRING "Nr. of processes used for testing")

set(TEST_OMP_THREADS "1" CACHE STRING "Nr. of OpeMP-threads used for testing")


# Installation paths
set(INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH
  "Installation directory for executables")

set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" CACHE PATH "Installation directory for libraries")

set(INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include/dftb+" CACHE PATH
  "Installation directory for header and include files")

set(INSTALL_MOD_DIR "${INSTALL_INC_DIR}/modfiles" CACHE PATH
  "Installation directory for Fortran module files")

option(EXPORT_EXTLIBS_WITH_PATH
  "Whether external libraries in the CMake export file should contain their full path" FALSE)
# For CMake experts only: It allows to link exact the same external libraries when using
# the library in an other CMake project. It does not play well with the CMake export file of
# the ELSI library.

set(PKGCONFIG_LANGUAGE "Fortran" CACHE STRING
  "Compiler and Linker language to assume when creating the pkg-config export file (C or Fortran)")
# The pkg-config export file (lib/pkgconfig/dftbplus.pc) contains the compiler and linker options
# needed to link the DFTB+ library to an application. (It can be queried with the pkg-config tool.)
# Depending on the language setting ("C" or "Fortran") you would get the flags for the case of using
# that compiler for the linking.


####################################################################################################
#
# NOTE FOR DEVELOPERS: Do not customise any settings here or in any of the sys/${ARCH}.cmake files
# as they contain the official defaults DFTB+ is shipped with. (Except you have a good reason to
# change such a default). If you need to customise any of the settings for your system, create a
# custom cmake file (e.g. custom.cmake) containing (only) the settings you would like to
# override. For an example, see
#
#     https://gist.github.com/aradi/39ab88acfbacc3b2f44d1e41e4da15e7
#
# When invoking CMake, pre-populate its cache with your custom settings using the -C option. For
# example, assuming the DFTB+ source is in ~/dftbplus, issue:
#
#     cmake -C ~/dftbplus/custom.cmake ~/dftbplus
#
# The settings in custom.cmake will take precedence over the corresponding settings in config.cmake
# and sys/${ARCH}.cmake. Make sure, you do *not* put your customised makefile under version control.
#
# Alternatively, you may also override settings on the command line, e.g.:
#
#     cmake -DWITH_SOCKETS=1 ~/dftbplus
#
####################################################################################################
