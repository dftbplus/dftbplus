cmake_minimum_required(VERSION 3.16)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include(DftbPlusUtils)

dftbp_ensure_out_of_source_build()
dftbp_load_build_settings()

set(DFTBPLUS_VERSION "22.1")

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

if(WITH_MPI)
  find_package(MPI REQUIRED)
  if(NOT MPI_FORTRAN_FOUND)
    message(FATAL_ERROR "Compiler ${CMAKE_Fortran_COMPILER} is not MPI capable but is specified "
      "for a WITH_MPI=TRUE build")
  endif()
endif()

if(WITH_OMP)
  if(NOT TARGET OpenMP::OpenMP_Fortran)
    find_package(OpenMP REQUIRED)
    # Fix CMake bug on OpenMP settings for the NAG compiler
    # See https://gitlab.kitware.com/cmake/cmake/-/issues/21280
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG" AND WITH_OMP)
      set_property(TARGET OpenMP::OpenMP_Fortran PROPERTY INTERFACE_LINK_LIBRARIES "")
      set_property(TARGET OpenMP::OpenMP_Fortran PROPERTY INTERFACE_LINK_OPTIONS "-openmp")
    endif()
  endif()
else()
  # Create empty targets, so that subprojects can depend on it without actual effect
  if(NOT TARGET OpenMP::OpenMP_Fortran)
    add_library(OpenMP::OpenMP_Fortran INTERFACE IMPORTED)
  endif()
  if(NOT TARGET OpenMP::OpenMP_C)
    add_library(OpenMP::OpenMP_C INTERFACE IMPORTED)
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

if(WITH_GPU AND NOT WITH_MPI)
  set(WITH_MAGMA TRUE)
else()
  set(WITH_MAGMA FALSE)
endif()

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

  if(WITH_GPU)
    find_package(CUDAToolkit REQUIRED)
  endif()

else()
  list(APPEND FYPP_FLAGS -DELSI_VERSION=0)
endif()

if(WITH_PLUMED)
  find_package(CustomPlumed REQUIRED)
  list(APPEND PKG_CONFIG_REQUIRES plumed)
endif()

find_package(CustomBlas REQUIRED)
find_package(CustomLapack REQUIRED)
list(APPEND PKG_CONFIG_LIBS_PRIVATE ${BLAS_LIBRARY} ${LAPACK_LIBRARY})

if(WITH_ARPACK)
  find_package(CustomArpack REQUIRED)
  list(APPEND PKG_CONFIG_LIBS_PRIVATE ${ARPACK_LIBRARY})
endif()

if(WITH_MAGMA)
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

# Include ddCOSMO library
add_subdirectory(external/ddcosmo EXCLUDE_FROM_ALL)

if(WITH_SOCKETS)
  add_subdirectory(external/fsockets EXCLUDE_FROM_ALL)
endif()

if(WITH_CHIMES)
  enable_language(CXX)
  set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 0)
endif()

# If INCLUDE_INDIRECT_DEPS is non-empty, indirect dependencies must also be explicitely treated
string(REGEX MATCH "(^|;)[Ss]ubmodule(^|;)" INCLUDE_INDIRECT_DEPS "${HYBRID_CONFIG_METHODS}")

# Note: GIT_TAG hashes below must be updated with the utils/test/check_submodule_commits script!

if(WITH_MPI)
  set(MPIFX_GIT_REPOSITORY "https://github.com/dftbplus/mpifx.git")
  set(MPIFX_GIT_TAG "da51073aa87831a91ae756a9773e39ea000d1c3a")  # do not change manually!
  dftbp_config_hybrid_dependency(MpiFx MpiFx::MpiFx "${HYBRID_CONFIG_METHODS}" "QUIET"
    external/mpifx "${exclude}" "${MPIFX_GIT_REPOSITORY}" "${MPIFX_GIT_TAG}")

  # Find ScaLAPACK to make sure, all components (e.g. ScalapackFx, MBD) use the same ScaLAPACK
  find_package(CustomScalapack REQUIRED)

  set(SCALAPACKFX_GIT_REPOSITORY "https://github.com/dftbplus/scalapackfx.git")
  set(SCALAPACKFX_GIT_TAG "7b98ba99d595a0760da12895b1a63959e881182c")  # do not change manually!
  dftbp_config_hybrid_dependency(ScalapackFx ScalapackFx::ScalapackFx "${HYBRID_CONFIG_METHODS}"
    "QUIET" external/scalapackfx "${exclude}" "${SCALAPACKFX_GIT_REPOSITORY}"
    "${SCALAPACKFX_GIT_TAG}")
  list(APPEND PKG_CONFIG_REQUIRES mpifx scalapackfx)
endif()

# MBD must be invoked *after* Scalapack has been found by ScalapackFx
if(WITH_MBD)
  set(MBD_GIT_REPOSITORY "https://github.com/libmbd/libmbd.git")
  set(MBD_GIT_TAG "95c95b5c1606b4576d49082a0396ea5701333375")  # do not change manually!
  dftbp_config_hybrid_dependency(Mbd Mbd::Mbd "${HYBRID_CONFIG_METHODS}" "QUIET" external/mbd
    "${exclude}" "${MBD_GIT_REPOSITORY}" "${MBD_GIT_TAG}")
  #list(APPEND PKG_CONFIG_REQUIRES mbd)
endif()

if(WITH_TRANSPORT)
  set(LIBNEGF_GIT_REPOSITORY "https://github.com/libnegf/libnegf.git")
  set(LIBNEGF_GIT_TAG "5b769f76406b12ec2faad30593cb061c5bec7d47")  # do not change manually!
  dftbp_config_hybrid_dependency(Negf Negf::Negf "${HYBRID_CONFIG_METHODS}" "QUIET"
    external/libnegf "${exclude}" "${LIBNEGF_GIT_REPOSITORY}" "${LIBNEGF_GIT_TAG}")
  #list(APPEND PKG_CONFIG_REQUIRES negf)
endif()
add_subdirectory(external/mudpack)

if(WITH_TBLITE OR WITH_SDFTD3)
  set(MCTC_LIB_GIT_REPOSITORY "https://github.com/grimme-lab/mctc-lib.git")
  set(MCTC_LIB_GIT_TAG "3d0e921085eb6ca535dfceec876bc99baf52f05d")  # do not change manually!
  dftbp_config_hybrid_dependency(mctc-lib mctc-lib::mctc-lib "${HYBRID_CONFIG_METHODS}" "QUIET"
    external/mctc-lib "${exclude}" "${MCTC_LIB_GIT_REPOSITORY}" "${MCTC_LIB_GIT_TAG}")
  #list(APPEND PKG_CONFIG_REQUIRES mctc-lib)

  if(INCLUDE_INDIRECT_DEPS)
    set(MSTORE_GIT_REPOSITORY "https://github.com/grimme-lab/mstore.git")
    set(MSTORE_GIT_TAG "90fc6fc4d8a6c3d72e5a910347a45edabf65343c")  # do not change manually!
    dftbp_config_hybrid_dependency(mstore mstore::mstore "${HYBRID_CONFIG_METHODS}" "QUIET"
      external/mstore "${exclude}" "${MSTORE_GIT_REPOSITORY}" "${MSTORE_GIT_TAG}")
    #list(APPEND PKG_CONFIG_REQUIRES mstore)
  endif()

  if(INCLUDE_INDIRECT_DEPS)
    set(TOML_F_GIT_REPOSITORY "https://github.com/toml-f/toml-f.git")
    set(TOML_F_GIT_TAG "f066ec6e7fb96d8faf83ab6614ee664a26ad8d57")  # do not change manually!
    dftbp_config_hybrid_dependency(toml-f toml-f::toml-f "${HYBRID_CONFIG_METHODS}" "QUIET"
      external/toml-f "${exclude}" "${TOML_F_GIT_REPOSITORY}" "${TOML_F_GIT_TAG}")
    #list(APPEND PKG_CONFIG_REQUIRES toml-f)
  endif()

  set(S_DFTD3_GIT_REPOSITORY "https://github.com/awvwgk/simple-dftd3.git")
  set(S_DFTD3_GIT_TAG "02dead293e59408100eb382e97be56e2a03fb317")  # do not change manually!
  dftbp_config_hybrid_dependency(s-dftd3 s-dftd3::s-dftd3 "${HYBRID_CONFIG_METHODS}" "QUIET"
    external/s-dftd3 "${exclude}" "${S_DFTD3_GIT_REPOSITORY}" "${S_DFTD3_GIT_TAG}")
  #list(APPEND PKG_CONFIG_REQUIRES s-dftd3)
endif()

if(WITH_TBLITE)
  if(INCLUDE_INDIRECT_DEPS)
    set(MULTICHARGE_GIT_REPOSITORY "https://github.com/grimme-lab/multicharge.git")
    set(MULTICHARGE_GIT_TAG "46ea4e370b1a59e01cc022d1f65c64bae8dbf968")  # do not change manually!
    dftbp_config_hybrid_dependency(multicharge multicharge::multicharge "${HYBRID_CONFIG_METHODS}"
      "QUIET" external/multicharge "${exclude}" "${MULTICHARGE_GIT_REPOSITORY}"
      "${MULTICHARGE_GIT_TAG}")
    #list(APPEND PKG_CONFIG_REQUIRES multicharge)
  endif()

  if(INCLUDE_INDIRECT_DEPS)
    set(DFTD4_GIT_REPOSITORY "https://github.com/dftd4/dftd4.git")
    set(DFTD4_GIT_TAG "f572723f5a5cb8243ea88b0253e182088311c6f3")  # do not change manually!
    dftbp_config_hybrid_dependency(dftd4 dftd4::dftd4 "${HYBRID_CONFIG_METHODS}" "QUIET"
      external/dftd4 "${exclude}" "${DFTD4_GIT_REPOSITORY}" "${DFTD4_GIT_TAG}")
    #list(APPEND PKG_CONFIG_REQUIRES dftd4)
  endif()

  set(TBLITE_GIT_REPOSITORY "https://github.com/tblite/tblite.git")
  set(TBLITE_GIT_TAG "0542ce7ae0e323941156949a0620ca260bc0ce7f")  # do not change manually!
  dftbp_config_hybrid_dependency(tblite tblite::tblite "${HYBRID_CONFIG_METHODS}" "QUIET"
    external/tblite "${exclude}" "${TBLITE_GIT_REPOSITORY}" "${TBLITE_GIT_TAG}")
  #list(APPEND PKG_CONFIG_REQUIRES tblite)

endif()

# ChIMES currently can only be integrated as a submodule
if(WITH_CHIMES)
  set(CHIMES_GIT_REPOSITORY "https://github.com/dftbplus/chimes_calculator.git")
  set(CHIMES_GIT_TAG "98e365365d3b3e4b11b363b4a50dd2ad6cc825d4")  # do not change manually!
  dftbp_config_hybrid_dependency(ChimesCalc ChimesCalc::ChimesCalc_Fortran
    "${HYBRID_CONFIG_METHODS}" "QUIET"
  external/chimes "${exclude}" "${CHIMES_GIT_REPOSITORY}" "${CHIMES_GIT_TAG}")
endif()

# FyTest currently can only be integrated as a submodule
if(NOT BUILD_EXPORTED_TARGETS_ONLY AND WITH_UNIT_TESTS)
  set(FYTEST_GIT_REPOSITORY "https://github.com/aradi/fytest.git")
  set(FYTEST_GIT_TAG "94e05acb31e78da1180845e345e542a1a28577d4")  # do not change manually!
  dftbp_config_hybrid_dependency(Fytest FyTest::FyTest "Submodule" "QUIET"
  external/fytest "${exclude}" "${FYTEST_GIT_REPOSITORY}" "${FYTEST_GIT_TAG}")
endif()


#
# Add internal components
#
add_subdirectory(src/dftbp)
add_subdirectory(app)
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


#
# Coverage testing related targets
#
if(LCOV_REPORT)
  find_program(LCOV lcov REQUIRED)
  find_program(GENHTML genhtml)
  dftbp_create_lcov_targets("${LCOV}" "${GENHTML}" "${CMAKE_CURRENT_BINARY_DIR}/lcov"
    "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/prog")
endif()
