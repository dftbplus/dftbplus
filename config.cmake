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

option(WITH_PEXSI "Whether the ELSI libraries are compiled with PEXSI support" FALSE)
# Works only with MPI-parallel build and ELSI enabled.

option(WITH_GPU "Whether DFTB+ should support GPU-acceleration via the MAGMA-library" FALSE)

option(WITH_TRANSPORT "Whether transport via libNEGF should be included." FALSE)

option(WITH_SOCKETS "Whether socket communication should be allowed for" FALSE)

option(WITH_ARPACK "Whether the ARPACK library should be included (needed for TD-DFTB)" FALSE)
# Works only with non-MPI (serial) build

option(WITH_DFTD3 "Whether the DFTD3 library should be included" FALSE)
# NOTE: Due to the license of the DFTD3 library, the combined code must be distributed under the
# GPLv3 license (as opposed to the LGPLv3 license of the DFTB+ package)

option(BUILD_API "Whether the high-level API to the DFTB+ library should be built" FALSE)
# Turn this on, if you want to use DFTB+ as a library (instead of a standalone program)

option(MONOLITHIC_LIBDFTBPLUS
  "Whether the DFTB+ library built should contain some of the libraries it depends on" FALSE)

#
# Architecture dependent build settings
#
set(ARCH "x86_64-linux-gnu" CACHE STRING
  "Selects which architecture dependent settings should be used")

# Include architecture dependant build settings from the sys-directory
include(${CMAKE_SOURCE_DIR}/sys/${ARCH}.cmake)


set(TEST_MPI_PROCS "1" CACHE STRING "Nr. of processes used for testing")

set(TEST_OMP_THREADS "1" CACHE STRING "Nr. of OpeMP-threads used for testing")


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
