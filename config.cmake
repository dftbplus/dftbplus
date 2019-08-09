#
# Global architecture independent build settings
#
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Release|Debug)")

set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/_install" CACHE STRING
  "Directory to install the compiled code into")

option(WITH_SOCKETS "Whether socket communication should be allowed for" FALSE)

option(WITH_DFTD3 "Whether the DFTD3 library should be included" FALSE)

option(WITH_ARPACK "Whether the ARPACK library should be included (needed for TD-DFTB)" FALSE)

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
# change such a default for some reasons). If you need to customise any of the settings for your
# system, create a custom cmake file (e.g. custom.cmake) containing (only) the settings you would
# like to override. When invoking CMake, pre-populate its cache with your custom settings using the
# -C option. For example, assuming the DFTB+ source is in ~/dftbplus, issue:
#
# cmake -C ~/dftbplus/custom.cmake ~/dftbplus
#
# The settings in custom.cmake will take precedence over the corresponding settings in config.cmake
# and sys/${ARCH}.cmake. Make sure, you do *not* put your customised makefile under version control.
#
# Alternatively, you may also override settings on the command line, e.g.:
#
# cmake -DWITH_SOCKETS=1 ~/dftbplus
#
####################################################################################################
