option(WITH_SOCKETS "Whether socket communication should allowed for" FALSE)

option(WITH_DFTD3 "Whether the DFTD3 library should be included" FALSE)

option(WITH_ARPACK "Whether the ARPACK library should be included" FALSE)

set(TEST_MPI_PROCS "1" CACHE STRING "Nr. of processes used for testing")

set(TEST_OMP_THREADS "1" CACHE STRING "Nr. of OpeMP-threads used for testing")

option(MONOLITHIC_LIBDFTBPLUS
  "Whether the DFTB+ library built should contain some of the libraries it depends on" FALSE)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Release|Debug)")

# Include optional architecture dependent manual settings
include(${CMAKE_SOURCE_DIR}/arch.cmake OPTIONAL)

