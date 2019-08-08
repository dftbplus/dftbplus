cmake_minimum_required(VERSION 3.5)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(dftbpCommon)

include(${CMAKE_SOURCE_DIR}/config.cmake)

project(dftbplus VERSION 19.2 LANGUAGES Fortran C)

set(FYPP "${PROJECT_SOURCE_DIR}/external/fypp/bin/fypp" CACHE FILEPATH "Fypp preprocessor")

set(PYTHON_INTERPRETER "python" CACHE STRING
  "Python interpreter to use for installing and test python components")

if(WITH_MPI)
  set(TEST_RUNNER "env OMP_NUM_THREADS=${TEST_OMP_THREADS} mpiexec -n ${TEST_MPI_PROCS}")
else()
  set(TEST_RUNNER "env OMP_NUM_THREADS=${TEST_OMP_THREADS}")
endif()


if(WITH_MPI)
  find_package(MPI REQUIRED)
endif()

find_package(OpenMP REQUIRED)
string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
string(APPEND CMAKE_Fortran_FLAGS " ${OpenMP_Fortran_FLAGS}")
string(APPEND CMAKE_EXE_LINKER_FLAGS " ${OpenMP_EXE_LINKER_FLAGS}")

#
# Check external libraries
#
if(WITH_MPI)
  create_library_targets("${SCALAPACK_LIBRARIES}" "${SCALAPACK_LIBRARY_PATHS}")
endif()

create_library_targets("${LAPACK_LIBRARIES}" "${LAPACK_LIBRARY_PATHS}")
create_library_targets("${OTHER_LIBRARIES}" "${OTHER_LIBRARY_PATHS}")

if(WITH_ARPACK)
  create_library_targets("${ARPACK_LIBRARIES}" "${ARPACK_LIBRARY_PATHS}")
endif()


set(INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include/dftb+")
set(INSTALL_MOD_DIR "${CMAKE_INSTALL_PREFIX}/include/dftb+/mods")

dftbp_get_release_name(RELEASE)

#
# Add optional external components
#
add_subdirectory(external/xmlf90 EXCLUDE_FROM_ALL)

if(WITH_SOCKETS)
  add_subdirectory(external/fsockets EXCLUDE_FROM_ALL)
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
