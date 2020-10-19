cmake_minimum_required(VERSION 3.16)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include(DftbPlusUtils)

dftbp_ensure_out_of_source_build()
dftbp_load_build_settings()

set(DFTBPLUS_VERSION "20.1")

project(dftbplus VERSION ${DFTBPLUS_VERSION} LANGUAGES Fortran C)

dftbp_setup_build_type()
dftbp_load_toolchain_settings()
dftbp_setup_global_compiler_flags()

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

find_package(MPI QUIET)
if(WITH_MPI AND NOT MPI_FORTRAN_FOUND)
  message(FATAL_ERROR "Compiler ${CMAKE_Fortran_COMPILER} is not MPI capable but is specified for "
    "a WITH_MPI=TRUE build")
elseif(NOT WITH_MPI AND MPI_FORTRAN_FOUND)
  message(WARNING "MPI enabled compiler ${CMAKE_Fortran_COMPILER} found for a non-MPI build. Your "
    "build will NOT be MPI-parallelised. Set WITH_MPI=TRUE in order to obtain an MPI-parallelised "
    "build.")
endif()

if(WITH_OMP AND NOT TARGET OpenMP::OpenMP_Fortran)
  find_package(OpenMP REQUIRED)
  # Fix CMake bug on OpenMP settings for the NAG compiler
  # See https://gitlab.kitware.com/cmake/cmake/-/issues/21280
  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG" AND WITH_OMP)
    set_property(TARGET OpenMP::OpenMP_Fortran PROPERTY INTERFACE_LINK_LIBRARIES "")
    set_property(TARGET OpenMP::OpenMP_Fortran PROPERTY INTERFACE_LINK_OPTIONS "-openmp")
  endif()
endif()

#
# Check external libraries
#
set(PKG_CONFIG_REQUIRES)
set(PKG_CONFIG_LIBS)
set(PKG_CONFIG_LIBS_PRIVATE)

if(WITH_ELSI AND NOT TARGET elsi::elsi)
  set(ELSI_MIN_VERSION "2.5.0")
  find_package(elsi ${ELSI_MIN_VERSION} REQUIRED PATHS ${ELSI_ROOT})
  list(APPEND PKG_CONFIG_REQUIRES elsi)

  if(TARGET elsi::pexsi)
    set(ELSI_WITH_PEXSI TRUE)
    enable_language(CXX)
    # Avoid to link everything associated with PEXSI with C++
    set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 0)
  else()
    set(ELSI_WITH_PEXSI FALSE)
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
  find_package(CustomPlumed REQUIRED)
  list(APPEND PKG_CONFIG_REQUIRES plumed)
endif()

find_package(CustomLapack REQUIRED)
list(APPEND PKG_CONFIG_LIBS_PRIVATE ${LAPACK_LIBRARY})

if(WITH_ARPACK)
  find_package(CustomArpack REQUIRED)
  list(APPEND PKG_CONFIG_LIBS_PRIVATE ${ARPACK_LIBRARY})
endif()

if(WITH_GPU)
  find_package(CustomMagma REQUIRED)
  list(APPEND PKG_CONFIG_REQUIRES magma)
endif()


#
# Preprocessor details
#
set(FYPP "${PROJECT_SOURCE_DIR}/external/fypp/bin/fypp" CACHE FILEPATH "Fypp preprocessor")
dftbp_add_fypp_defines(FYPP_FLAGS)

set(FYPP_CONFIG_FLAGS "${FYPP_FLAGS}")
# Make sure, the line-marker option is not set
list(REMOVE_ITEM FYPP_CONFIG_FLAGS "-n")
set(FYPP_BUILD_FLAGS "${FYPP_FLAGS}" "$<IF:$<CONFIG:Debug>,-DDEBUG=1,-DDEBUG=0>")

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

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

add_subdirectory(external/xmlf90 EXCLUDE_FROM_ALL)

if(WITH_SOCKETS)
  add_subdirectory(external/fsockets EXCLUDE_FROM_ALL)
endif()

if(WITH_DFTD3)
  add_subdirectory(external/dftd3 ${exclude})
endif()

if(WITH_MPI)
  dftbp_build_or_find_external(MpiFx BUILD_MPIFX external/mpifx "${exclude}")
  dftbp_build_or_find_external(ScalapackFx BUILD_SCALAPACKFX external/scalapackfx "${exclude}")
endif()

# MBD must be invoked *after* Scalapack has been found by ScalapackFx
if(WITH_MBD)
  dftbp_build_or_find_external(Mbd BUILD_MBD external/mbd "${exclude}")
endif()

if(WITH_TRANSPORT)
  dftbp_build_or_find_external(Negf BUILD_NEGF external/libnegf "${exclude}")
endif()
add_subdirectory(external/mudpack)


#
# Add internal components
#
add_subdirectory(prog)
#add_subdirectory(tools)


#
# Testing
#
string(CONFIGURE "${TEST_RUNNER_TEMPLATE}" TEST_RUNNER)
if(NOT WITH_MPI)
  string(CONFIGURE "${MODES_RUNNER_TEMPLATE}" MODES_RUNNER)
endif()
enable_testing()
if(NOT BUILD_EXPORTED_TARGETS_ONLY)
  add_subdirectory(test)
endif()


#
# Export package info
# 
if(WITH_API)
  include(CMakePackageConfigHelpers)

  add_library(DftbPlus INTERFACE)
  target_link_libraries(DftbPlus INTERFACE dftbplus)
  install(TARGETS DftbPlus EXPORT dftbplus-targets)

  install(EXPORT dftbplus-targets
    FILE dftbplus-targets.cmake
    NAMESPACE DftbPlus::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/dftbplus")
  
  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/export/dftbplus-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/dftbplus-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dftbplus)

  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/dftbplus-config-version.cmake
    VERSION ${API_VERSION}
    COMPATIBILITY SameMajorVersion)

  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/dftbplus-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/dftbplus-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dftbplus)

  install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dftbplus)

  dftbp_get_pkgconfig_params(PKGCONFIG_REQUIRES PKGCONFIG_LIBS PKGCONFIG_LIBS_PRIVATE
    PKGCONFIG_C_FLAGS)
  configure_file(utils/export/dftbplus.pc.in ${CMAKE_CURRENT_BINARY_DIR}/dftbplus.pc @ONLY)
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/dftbplus.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
  
endif()
