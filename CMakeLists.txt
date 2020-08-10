cmake_minimum_required(VERSION 3.5)
cmake_policy(VERSION 3.5)

set(DFTBPLUS_VERSION "20.1")

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(dftbpUtils)

dftbp_ensure_out_of_source_build()
dftbp_load_build_settings()

project(dftbplus VERSION ${DFTBPLUS_VERSION} LANGUAGES Fortran C)

# By default, use intrinsic Fortran 2008 erf/erfc
set(INTERNAL_ERFC CACHE BOOL 0)

dftbp_load_toolchain_settings()
dftbp_setup_global_compiler_flags()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Neither CMake nor we have any ideas, how gcc could link NAG-compiled Fortran code with OMP.
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG" AND WITH_OMP)
  set(WITH_C_EXECUTABLES FALSE)
else()
  set(WITH_C_EXECUTABLES TRUE)
endif()

dftbp_ensure_config_consistency()

dftbp_get_release_name(RELEASE)

if(WITH_API)
  dftbp_get_api_version(API_VERSION API_VERSION_MAJOR API_VERSION_MINOR API_VERSION_PATCH)
endif()

find_package(MPI)
if(WITH_MPI AND NOT MPI_FORTRAN_FOUND)
  message(FATAL_ERROR "Compiler ${CMAKE_Fortran_COMPILER} is not MPI capable but is specified for "
    "a WITH_MPI=TRUE build")
elseif(NOT WITH_MPI AND MPI_FORTRAN_FOUND)
  message(WARNING "MPI enabled compiler ${CMAKE_Fortran_COMPILER} found for a non-MPI build. Your "
    "build will NOT be MPI-parallelised. Set WITH_MPI=TRUE in order to obtain an MPI-parallelised "
    "build.")
endif()

if(WITH_OMP)
  find_package(OpenMP REQUIRED)
  string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
  string(APPEND CMAKE_Fortran_FLAGS " ${OpenMP_Fortran_FLAGS}")
  string(APPEND CMAKE_EXE_LINKER_FLAGS " ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#
# Check external libraries
#
set(EXPORTED_EXTERNAL_LIBRARIES)
set(EXPORTED_EXTERNAL_LIBRARY_DIRS)
set(EXPORTED_EXTERNAL_PACKAGES)
set(EXPORTED_EXPORTED_COMPILED_LIBRARIES)

if(WITH_MPI)
  dftbp_create_library_targets("${SCALAPACK_LIBRARIES}" "${SCALAPACK_LIBRARY_DIRS}")
  list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${SCALAPACK_LIBRARIES})
  list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${SCALAPACK_LIBRARY_DIRS})
endif()

if(WITH_ELSI)
  set(ELSI_MIN_VERSION "2.5.0")
  find_package(elsi ${ELSI_MIN_VERSION} REQUIRED PATHS ${ELSI_ROOT})
  dftbp_create_library_targets("${ELSI_EXTERNAL_LIBRARIES}" "${ELSI_EXTERNAL_LIBRARY_DIRS}")
  list(APPEND EXPORTED_EXTERNAL_PACKAGES elsi)
  list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${ELSI_EXTERNAL_LIBRARIES})
  list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${ELSI_EXTERNAL_LIBRARY_DIRS})
  if(TARGET elsi::pexsi)
    set(ELSI_WITH_PEXSI TRUE)
  else()
    set(ELSI_WITH_PEXSI FALSE)
  endif()
  if(ELSI_WITH_PEXSI)
    dftbp_create_library_targets("${PEXSI_EXTERNAL_LIBRARIES}" "${PEXSI_EXTERNAL_LIBRARY_DIRS}")
    list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${PEXSI_EXTERNAL_LIBRARIES})
    list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${PEXSI_EXTERNAL_LIBRARY_DIRS})
  endif()
  if(elsi_VERSION VERSION_GREATER "2.5.0")
    message(STATUS "ELSI > v2.5.0")
    list(APPEND FYPP_FLAGS -DELSI_VERSION=2.6)
  else()
    list(APPEND FYPP_FLAGS -DELSI_VERSION=2.5)
    message(STATUS "ELSI v2.5.0")
  endif()
else()
  list(APPEND FYPP_FLAGS -DELSI_VERSION=0)
endif()

if(WITH_PLUMED)
  #find_package(PkgConfig REQUIRED)
  #pkg_check_modules(PLUMED REQUIRED IMPORTED_TARGET plumed)
  dftbp_create_library_targets("${PLUMED_LIBRARIES}" "${PLUMED_LIBRARY_DIRS}")
  list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${PLUMED_LIBRARIES})
  list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${PLUMED_LIBRARY_DIRS})
endif()

dftbp_create_library_targets("${LAPACK_LIBRARIES}" "${LAPACK_LIBRARY_DIRS}")
list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${LAPACK_LIBRARIES})
list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${LAPACK_LIBRARY_DIRS})

if(WITH_ARPACK)
  dftbp_create_library_targets("${ARPACK_LIBRARIES}" "${ARPACK_LIBRARY_DIRS}")
  list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${ARPACK_LIBRARIES})
  list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${ARPACK_LIBRARY_DIRS})
endif()

if(WITH_GPU)
  find_package(MAGMA REQUIRED)
  dftbp_create_library_targets("${MAGMA_LIBRARIES}" "${MAGMA_LIBRARY_DIRS}")
  list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${MAGMA_LIBRARIES})
  list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${MAGMA_LIBRARY_DIRS})
endif()

dftbp_create_library_targets("${OTHER_LIBRARIES}" "${OTHER_LIBRARY_DIRS}")
list(APPEND EXPORTED_EXTERNAL_LIBRARIES ${OTHER_LIBRARIES})
list(APPEND EXPORTED_EXTERNAL_LIBRARY_DIRS ${OTHER_LIBRARY_DIRS})

#
# Preprocessor details
#
set(FYPP "${PROJECT_SOURCE_DIR}/external/fypp/bin/fypp" CACHE FILEPATH "Fypp preprocessor")
dftbp_add_fypp_defines(FYPP_FLAGS)

set(PYTHON_INTERPRETER "python3" CACHE STRING
  "Python interpreter to use for installing and test python components")


#
# Add optional external components
#
if(WITH_API)
  set(exclude)
  option(INSTALL_INCLUDE_FILES "Whether include files should be installed" TRUE)
elseif(BUILD_SHARED_LIBS)
  set(exclude)
  option(INSTALL_INCLUDE_FILES "Whether include files should be installed" FALSE)
else()
  set(exclude EXCLUDE_FROM_ALL)
  option(INSTALL_INCLUDE_FILES "Whether include files should be installed" FALSE)
endif()

add_subdirectory(external/xmlf90 EXCLUDE_FROM_ALL)

if(WITH_SOCKETS)
  add_subdirectory(external/fsockets EXCLUDE_FROM_ALL)
endif()

if(WITH_DFTD3)
  add_subdirectory(external/dftd3 ${exclude})
endif()

if(WITH_MBD)
  add_subdirectory(external/mbd ${exclude})
endif()

if(WITH_MPI)
  add_subdirectory(external/mpifx ${exclude})
  add_subdirectory(external/scalapackfx ${exclude})
endif()

if(WITH_TRANSPORT)
  add_subdirectory(external/libnegf ${exclude})
endif()
add_subdirectory(external/mudpack)

#
# Add internal components
#
add_subdirectory(prog)
add_subdirectory(tools)

#
# Testing
#
string(CONFIGURE "${TEST_RUNNER_TEMPLATE}" TEST_RUNNER)
if(NOT WITH_MPI)
  string(CONFIGURE "${MODES_RUNNER_TEMPLATE}" MODES_RUNNER)
endif()
enable_testing()
add_subdirectory(test)


#
# Export package info
# 
if(WITH_API)

  install(EXPORT dftbplus-targets FILE DftbPlusTargets.cmake NAMESPACE DftbPlus::
    DESTINATION ${INSTALL_LIB_DIR}/cmake/DftbPlus)
  
  configure_file(utils/export/DftbPlusConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/DftbPlusConfig.cmake @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/DftbPlusConfig.cmake DESTINATION
    ${INSTALL_LIB_DIR}/cmake/DftbPlus)

  dftbp_get_pkgconfig_params(PKGCONFIG_REQUIRES PKGCONFIG_LIBS PKGCONFIG_LIBS_PRIVATE
    PKGCONFIG_C_FLAGS)
  configure_file(utils/export/dftbplus.pc.in ${CMAKE_CURRENT_BINARY_DIR}/dftbplus.pc @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/dftbplus.pc DESTINATION ${INSTALL_LIB_DIR}/pkgconfig)
  
endif()
