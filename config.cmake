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

option(WITH_SOCKETS "Whether socket communication should be allowed for" FALSE)

option(WITH_ARPACK "Whether the ARPACK library should be included (needed for TD-DFTB)" FALSE)
# Works only with non-MPI (serial) build

option(WITH_DFTD3 "Whether the DFTD3 library should be included" FALSE)
# NOTE: Due to the license of the DFTD3 library, the combined code must be distributed under the
# GPLv3 license (as opposed to the LGPLv3 license of the DFTB+ package)

option(BUILD_API "Whether DFTB+ library with high-level API should be built and installed" FALSE)
# Turn this on, if you want to use DFTB+ as a library (libdftbplus.a) in order to integrate it into
# other software packages. (Otherwise only a stripped down version of the library without the public
# API would be built, and neither the library nor the module files would be installed.) When
# enabled, also the external libraries will be installed, which were compiled during the build
# process and are need when linking DFTB+ with other applications.


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
