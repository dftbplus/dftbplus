cmake_minimum_required(VERSION 3.16)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include(DftbPlusUtils)

dftbp_ensure_out_of_source_build()
dftbp_load_build_settings()

set(DFTBPLUS_VERSION "20.2.1")

project(dftbplus VERSION ${DFTBPLUS_VERSION} LANGUAGES Fortran C)

dftbp_setup_build_type()
dftbp_load_toolchain_settings()
dftbp_setup_global_compiler_flags()

# CMake has difficulties linking C-executables against NAG-compiled Fortran libraries
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
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
# Git repository check
#
find_package(Git)
if(Git_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  set(GIT_WORKING_COPY True)
endif()

#
# Check external libraries
#
set(PKG_CONFIG_REQUIRES)
set(PKG_CONFIG_LIBS)
set(PKG_CONFIG_LIBS_PRIVATE)

if(WITH_ELSI AND NOT TARGET elsi::elsi)
  set(ELSI_MIN_VERSION "2.5.0")
  find_package(elsi ${ELSI_MIN_VERSION} REQUIRED)
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

# Note: GIT_TAG hashes below must be updated with the utils/test/check_submodule_commits script!

if(WITH_MPI)
  set(MPIFX_GIT_REPOSITORY "https://github.com/dftbplus/mpifx.git")
  set(MPIFX_GIT_TAG "0cb07ee08cbb20f3f7bb2527152a4ec317c579ad")  # do not change manually!
  dftbp_config_hybrid_dependency(MpiFx MpiFx::MpiFx "${HYBRID_CONFIG_METHODS}" "QUIET" 
    external/mpifx "${exclude}" "${MPIFX_GIT_REPOSITORY}" "${MPIFX_GIT_TAG}")

  # Find ScaLAPACK to make sure, all components (e.g. ScalapackFx, MBD) use the same ScaLAPACK
  find_package(CustomScalapack REQUIRED)
  
  set(SCALAPACKFX_GIT_REPOSITORY "https://github.com/dftbplus/scalapackfx.git")
  set(SCALAPACKFX_GIT_TAG "cd7c2a789e7788fc42ea97cbfe92b5ea15f84ab7")  # do not change manually!
  dftbp_config_hybrid_dependency(ScalapackFx ScalapackFx::ScalapackFx "${HYBRID_CONFIG_METHODS}"
    "QUIET" external/scalapackfx "${exclude}" "${SCALAPACKFX_GIT_REPOSITORY}"
    "${SCALAPACKFX_GIT_TAG}")
  list(APPEND PKG_CONFIG_REQUIRES mpifx scalapackfx)
endif()

# MBD must be invoked *after* Scalapack has been found by ScalapackFx
if(WITH_MBD)
  set(MBD_GIT_REPOSITORY "https://github.com/jhrmnn/libmbd.git")
  set(MBD_GIT_TAG "82005cbb65bdf5d32ca021848eec8f19da956a77")  # do not change manually!
  dftbp_config_hybrid_dependency(Mbd Mbd::Mbd "${HYBRID_CONFIG_METHODS}" "QUIET" external/mbd
    "${exclude}" "${MBD_GIT_REPOSITORY}" "${MBD_GIT_TAG}")
  #list(APPEND PKG_CONFIG_REQUIRES mbd)
endif()

if(WITH_TRANSPORT)
  set(LIBNEGF_GIT_REPOSITORY "https://github.com/libnegf/libnegf.git")
  set(LIBNEGF_GIT_TAG "83d7588d3b05cf4aaf989085cea96f5dcf1b9830")  # do not change manually!
  dftbp_config_hybrid_dependency(Negf Negf::Negf "${HYBRID_CONFIG_METHODS}" "QUIET"
    external/libnegf "${exclude}" "${LIBNEGF_GIT_REPOSITORY}" "${LIBNEGF_GIT_TAG}")
  #list(APPEND PKG_CONFIG_REQUIRES negf)
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
