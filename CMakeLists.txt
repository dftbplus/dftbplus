cmake_minimum_required(VERSION 3.5)

project(dftbplus VERSION 19.2 LANGUAGES Fortran C)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(dftbpCommon)

set(SETTINGS_LOG "${CMAKE_BINARY_DIR}/dftbplus-settings.cmake")
dftbp_reset_file(${SETTINGS_LOG})

option(WITH_SOCKETS "Whether socket communication should allowed for" FALSE)
dftbp_log_settings(${SETTINGS_LOG} WITH_SOCKETS)

option(WITH_DFTD3 "Whether the DFTD3 library should be included" FALSE)
dftbp_log_settings(${SETTINGS_LOG} WITH_DFTD3)

option(WITH_ARPACK "Whether the ARPACK library should be included" FALSE)
dftbp_log_settings(${SETTINGS_LOG} WITH_ARPACK)

set(SCALAPACK_LIBRARIES "scalapack" CACHE STRING "ScaLAPACK library to link")
dftbp_log_settings(${SETTINGS_LOG} SCALAPACK_LIBRARIES)

set(FYPP "${PROJECT_SOURCE_DIR}/external/fypp/bin/fypp" CACHE FILEPATH "Fypp preprocessor")

set(TEST_MPI_PROCS "1" CACHE STRING "Nr. of processes used for testing")

set(TEST_OMP_THREADS "1" CACHE STRING "Nr. of OpeMP-threads used for testing")

set(PYTHON_INTERPRETER "python" CACHE STRING
  "Python interpreter to use for installing and test python components")

if(WITH_MPI)
  set(TEST_RUNNER "env OMP_NUM_THREADS=${TEST_OMP_THREADS} mpiexec -n ${TEST_MPI_PROCS}")
else()
  set(TEST_RUNNER "env OMP_NUM_THREADS=${TEST_OMP_THREADS}")
endif()
dftbp_log_settings(${SETTINGS_LOG} TEST_RUNNER)

find_package(LAPACK REQUIRED)
dftbp_log_settings(${SETTINGS_LOG} LAPACK_LIBRARIES)

if(WITH_MPI)
  find_package(MPI REQUIRED)
  #include_directories(SYSTEM ${MPI_Fortran_INCLUDE_PATH})
endif()

if(WITH_ARPACK)
  find_package(ARPACK REQUIRED)
  dftbp_log_settings(${SETTINGS_LOG} ARPACK_LIBRARIES)
endif()

find_package(OpenMP REQUIRED)
string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
string(APPEND CMAKE_Fortran_FLAGS " ${OpenMP_Fortran_FLAGS}")
string(APPEND CMAKE_EXE_LINKER_FLAGS " ${OpenMP_EXE_LINKER_FLAGS}")

set(INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include/dftb+")
set(INSTALL_MOD_DIR "${CMAKE_INSTALL_PREFIX}/include/dftb+/mods")

dftbp_get_release_name(RELEASE)


#
# Add optional external components
#
add_subdirectory(external/xmlf90)

if(WITH_SOCKETS)
  add_subdirectory(external/fsockets)
endif()

if(WITH_DFTD3)
  add_subdirectory(external/dftd3/origin EXCLUDE_FROM_ALL)
endif()

if(WITH_MPI)
  add_subdirectory(external/mpifx/origin EXCLUDE_FROM_ALL)
  add_subdirectory(external/scalapackfx/origin EXCLUDE_FROM_ALL)
endif()


#
# Add internal components
#
add_subdirectory(prog)
add_subdirectory(tools/dptools)

enable_testing()
add_subdirectory(test)
