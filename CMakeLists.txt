cmake_minimum_required(VERSION 3.5)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(dftbpUtils)

include(${CMAKE_SOURCE_DIR}/config.cmake)

project(dftbplus VERSION 19.2 LANGUAGES Fortran C)

dftbp_ensure_config_consistency()

set(FYPP "${PROJECT_SOURCE_DIR}/external/fypp/bin/fypp" CACHE FILEPATH "Fypp preprocessor")
dftbp_convert_to_list(FYPP_FLAGS)
dftbp_add_fypp_defines(FYPP_FLAGS)

set(PYTHON_INTERPRETER "python" CACHE STRING
  "Python interpreter to use for installing and test python components")

if(WITH_MPI)
  find_package(MPI REQUIRED)
endif()

if(WITH_OMP)
  find_package(OpenMP REQUIRED)
  string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
  string(APPEND CMAKE_Fortran_FLAGS " ${OpenMP_Fortran_FLAGS}")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " ${OpenMP_EXE_LINKER_FLAGS}")
endif()

dftbp_get_release_name(RELEASE)


#
# Check external libraries
#
dftbp_convert_to_list(SCALAPACK_LIBRARIES SCALAPACK_LIBRARY_DIRS LAPACK_LIBRARIES
  LAPACK_LIBRARY_DIRS ARPACK_LIBRARIES ARPACK_LIBRARY_DIRS ELSI_LIBRARIES ELSI_LIBRARY_DIRS
  ELSI_INCLUDE_DIRS PEXSI_LIBRARIES PEXSI_LIBRARY_DIRS PEXSI_INCLUDE_DIRS OTHER_LIBRARIES
  OTHER_LIBRARY_DIRS)

if(WITH_MPI)
  dftbp_create_library_targets("${SCALAPACK_LIBRARIES}" "${SCALAPACK_LIBRARY_DIRS}")
endif()

if(WITH_ELSI)
  dftbp_create_library_targets("${ELSI_LIBRARIES}" "${ELSI_LIBRARY_DIRS}")
endif()

if(WITH_PEXSI)
  dftbp_create_library_targets("${PEXSI_LIBRARIES}" "${PEXSI_LIBRARY_DIRS}")
endif()

dftbp_create_library_targets("${LAPACK_LIBRARIES}" "${LAPACK_LIBRARY_DIRS}")

if(WITH_ARPACK)
  dftbp_create_library_targets("${ARPACK_LIBRARIES}" "${ARPACK_LIBRARY_DIRS}")
endif()

dftbp_create_library_targets("${OTHER_LIBRARIES}" "${OTHER_LIBRARY_DIRS}")


#
# Add optional external components
#
add_subdirectory(external/xmlf90 EXCLUDE_FROM_ALL)

if(WITH_SOCKETS)
  add_subdirectory(external/fsockets EXCLUDE_FROM_ALL)
endif()

if(WITH_DFTD3)
  add_subdirectory(external/dftd3)
endif()

if(WITH_MPI)
  add_subdirectory(external/mpifx)
  add_subdirectory(external/scalapackfx)
endif()

if(WITH_TRANSPORT)
  if(BUILD_API)
    add_subdirectory(external/poisson)
    add_subdirectory(external/libnegf)
  else()
    add_subdirectory(external/poisson EXCLUDE_FROM_ALL)
    add_subdirectory(external/libnegf EXCLUDE_FROM_ALL)
  endif()
endif()

if(WITH_GPU)
  add_subdirectory(external/magmahelper EXCLUDE_FROM_ALL)
endif()


#
# Add internal components
#
add_subdirectory(prog)
add_subdirectory(tools/dptools)

enable_testing()
add_subdirectory(test)
