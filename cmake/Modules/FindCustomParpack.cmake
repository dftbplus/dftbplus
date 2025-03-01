# Distributed under the OSI-approved BSD 3-Clause License.
#
# Copyright (C) 2025  DFTB+ developers group
#

#[=======================================================================[.rst:
FindCustomParpack
-----------------

Finds the PARPACK library


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``PARPACK::PARPACK``
  The PARPACK library


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variables:

``CustomParpack_FOUND``
  True, if PARPACK had been found.

``PARPACK_FOUND``
  True if the system has the PARPACK library


Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``PARPACK_LIBRARY``
  Customized PARPACK library/libraries to use (instead of autodetected one). If
  not set or empty, the default library (parpack) will be tried.

``PARPACK_INCLUDE_DIRECTORY``
  Customized PARPACK include library to check for include files.

#]=======================================================================]

set(_parpack_found_method)
set(_parpack_found_msg_msg)

if (TARGET PARPACK::PARPACK)

  set(PARPACK_FOUND TRUE)
  set(_parpack_found_method "target")
  set(_parpack_found_msg "target PARPACK::PARPACK already defined")

elseif ("${CustomParpack_FOUND}" STREQUAL "")

  set(_quiet)
  if (CustomParpack_FIND_QUIETLY)
    set(_quiet QUIET)
  endif ()

  set(_required)
  if (CustomParpack_FIND_REQUIRED)
    set(_required REQUIRED)
  endif ()

  # Find dependency
  if (NOT TARGET BLAS::BLAS)
    find_package(BLAS ${_required} ${_quiet})
  endif ()

  # Use user specified value, if present
  if (PARPACK_LIBRARY)
    add_library(PARPACK::PARPACK UNKNOWN IMPORTED)
    set_target_properties(
      PARPACK::PARPACK PROPERTIES
      IMPORTED_LOCATION ${PARPACK_LIBRARY}
    )
    set(_parpack_found_method "user")
    set(_parpack_found_msg "as user defined library ${PARPACK_LIBRARY}")
    if (PARPACK_INCLUDE_DIRECTORY)
      set_target_properties(
        PARPACK::PARPACK PROPERTIES
        INCLUDE_DIRECTORIES ${PARPACK_INCLUDE_DIRECTORY}
      )
      set(_parpack_found_msg "${_parpack_found_msg} with user defined include dir ${PARPACK_INCLUDE_DIRECTORY}")
    endif ()
  endif ()

  # NOTE: Disabled as many ARPACK-NG installations contain broken CMake export files.
  # # Try to find it via CMake export file if had not been tried yet
  # if (NOT (arpackng_FOUND OR arpack-ng_FOUND))
  #   foreach (_package_name IN ITEMS arpackng arpack-ng)
  #     if (NOT TARGET PARPACK::PARPACK)
  #       find_package(${_package_name} QUIET)
  #       if (TARGET PARPACK::PARPACK)
  #         set(_parpack_found_method "cmake")
  #         set(_parpack_found_msg "as CMake package ${_package_name}")
  #       endif ()
  #     endif ()
  #     unset(_package_name)
  #   endforeach ()
  # endif ()

  # Try to find it via PkgConfig
  find_package(PkgConfig QUIET)
  if (PkgConfig_FOUND)
    if (NOT TARGET PARPACK::PARPACK)
      pkg_check_modules(parpack IMPORTED_TARGET parpack ${_quiet})
      if (TARGET PkgConfig::parpack)
        add_library(PARPACK::PARPACK INTERFACE IMPORTED)
        target_link_libraries(PARPACK::PARPACK INTERFACE PkgConfig::parpack)
        if (TARGET PARPACK::PARPACK)
          set(_parpack_found_method "pkgconfig")
          set(_parpack_found_msg "via package-config")
        endif ()
      endif ()
    endif ()
  endif ()

  # Try to find it by name
  if (NOT TARGET PARPACK::PARPACK)
    find_library(PARPACK_LIBRARIES parpack ${_quiet})
    find_path(PARPACK_INCLUDE_DIRECTORIES parpack.h HINTS PARPACK_INCLUDE_DIRECTORY ${_quiet})
    if (PARPACK_LIBRARIES)
      add_library(PARPACK::PARPACK INTERFACE IMPORTED)
      target_link_libraries(PARPACK::PARPACK INTERFACE ${PARPACK_LIBRARIES})
      set(_parpack_found_method "library")
      set(_parpack_found_msg "as library ${PARPACK_LIBRARIES}")
      if (PARPACK_INCLUDE_DIRECTORIES)
        set_target_properties(
          PARPACK::PARPACK PROPERTIES
          INCLUDE_DIRECTORIES ${PARPACK_INCLUDE_DIRECTORIES}
        )
        set(
          _parpack_found_msg
          "${_parpack_found_msg} with include directory ${PARPACK_INCLUDE_DIRECTORIES}"
        )
      endif ()
    endif ()
  endif ()

  if (TARGET PARPACK::PARPACK)
    set(PARPACK_FOUND TRUE)

    # Find and link dependency if not found via CMake export file
    if (NOT _parpack_found_method STREQUAL "cmake")
      if (NOT TARGET ARPACK::ARPACK)
        find_package(CustomArpack REQUIRED ${_quiet})
      endif ()
      target_link_libraries(PARPACK::PARPACK INTERFACE ARPACK::ARPACK)
    endif()
  endif()

  unset(_quiet)
  unset(_required)

endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CustomParpack DEFAULT_MSG _parpack_found_msg)

unset(_parpack_found_method)
unset(_parpack_found_msg)
