cmake_minimum_required(VERSION 3.5)

project(dftbplus VERSION 18.2 LANGUAGES Fortran C)

set(USER_CMAKE_DIR ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_MODULE_PATH ${USER_CMAKE_DIR})
include(${USER_CMAKE_DIR}/common.cmake)

set(SETTINGS_LOG "${CMAKE_BINARY_DIR}/dftbplus-settings.cmake")
_reset_file(${SETTINGS_LOG})

set(WITH_SOCKETS FALSE CACHE BOOL "Whether socket communication should allowed for")
_log_settings(${SETTINGS_LOG} WITH_SOCKETS)

set(WITH_DFTD3 FALSE CACHE BOOL "Whether the DFTD3 library should be included")
_log_settings(${SETTINGS_LOG} WITH_DFTD3)

set(WITH_ARPACK FALSE CACHE BOOL "Whether the ARPACK library should be included")
_log_settings(${SETTINGS_LOG} WITH_ARPACK)

set(FYPP "${PROJECT_SOURCE_DIR}/external/fypp/bin/fypp" CACHE FILEPATH "Fypp preprocessor")

find_package(LAPACK REQUIRED)
_log_settings(${SETTINGS_LOG} LAPACK_LIBRARIES)

if (WITH_ARPACK)
  find_package(ARPACK REQUIRED)
  _log_settings(${SETTINGS_LOG} ARPACK_LIBRARIES)
endif()

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include")
set(INSTALL_MOD_DIR "${CMAKE_INSTALL_PREFIX}/include")

_get_release_name(RELEASE)


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


#
# Add internal components
#
add_subdirectory(prog/dftb+)
add_subdirectory(prog/modes)
add_subdirectory(prog/waveplot)
