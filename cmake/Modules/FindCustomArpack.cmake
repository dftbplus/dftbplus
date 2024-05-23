# Distributed under the OSI-approved BSD 3-Clause License.
#
# Copyright (C) 2024  DFTB+ developers group
#

#[=======================================================================[.rst:
FindCustomArpack
----------------

Finds the ARPACK library


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``ARPACK::ARPACK``
  The ARPACK library


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variables:

``CustomArpack_FOUND``
  True, if ARPACK had been found.

``ARPACK_FOUND``
  True if the system has the ARPACK library


Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``ARPACK_LIBRARY``
  Customized ARPACK library/libraries to use (instead of autodetected one). If
  not set or empty, the default library (arpack) will be tried.

``ARPACK_INCLUDE_DIRECTORY``
  Customized ARPACK include library to check for include files.

#]=======================================================================]

set(_arpack_found_method)
set(_arpack_found_msg)

if (TARGET ARPACK::ARPACK)

  set(ARPACK_FOUND TRUE)
  set(_arpack_found_method "target")
  set(_arpack_found_msg "target ARPACK::ARPACK already defined")

elseif ("${CustomArpack_FOUND}" STREQUAL "")

  set(_quiet)
  if (CustomArpack_FIND_QUIETLY)
    set(_quiet QUIET)
  endif ()

  set(_required)
  if (CustomArpack_FIND_REQUIRED)
    set(_required REQUIRED)
  endif ()

  # Use user specified value, if present
  if (ARPACK_LIBRARY)
    add_library(ARPACK::ARPACK UNKNOWN IMPORTED)
    set_target_properties(
      ARPACK::ARPACK PROPERTIES
      IMPORTED_LOCATION ${ARPACK_LIBRARY}
    )
    set(_arpack_found_msg "as user defined library ${ARPACK_LIBRARY}")
    set(_arpack_found_method "user")
    if (ARPACK_INCLUDE_DIRECTORY)
      set_target_properties(
        ARPACK::ARPACK PROPERTIES
        INCLUDE_DIRECTORIES ${ARPACK_INCLUDE_DIRECTORY}
      )
      set(_arpack_found_msg "${_arpack_found_msg} with user defined include dir ${ARPACK_INCLUDE_DIRECTORY}")
    endif ()
  endif ()

  # NOTE: Disabled as many ARPACK-NG installations contain broken CMake export files.
  # # Try to find it via CMake export file if had not been tried yet
  # if (NOT (arpackng_FOUND OR arpack-ng_FOUND))
  #   foreach (_package_name IN ITEMS arpackng arpack-ng)
  #     if (NOT TARGET ARPACK::ARPACK)
  #       find_package(${_package_name} QUIET)
  #       if (TARGET ARPACK::ARPACK)
  #         set(_arpack_found_method "cmake")
  #         set(_arpack_found_msg "as CMake package ${_package_name}")
  #       endif ()
  #     endif ()
  #     unset(_package_name)
  #   endforeach ()
  # endif ()

  # Try to find it via PkgConfig
  find_package(PkgConfig QUIET)
  if (PkgConfig_FOUND)
    if (NOT TARGET ARPACK::ARPACK)
      pkg_check_modules(arpack IMPORTED_TARGET arpack ${_quiet})
      if (TARGET PkgConfig::arpack)
        add_library(ARPACK::ARPACK INTERFACE IMPORTED)
        target_link_libraries(ARPACK::ARPACK INTERFACE PkgConfig::arpack)
        if (TARGET ARPACK::ARPACK)
          set(_arpack_found_method "pkgconfig")
          set(_arpack_found_msg "via package-config")
        endif ()
      endif ()
    endif ()
  endif ()

  # Try to find it by name
  if (NOT TARGET ARPACK::ARPACK)
    find_library(ARPACK_LIBRARIES arpack ${_quiet})
    find_path(ARPACK_INCLUDE_DIRECTORIES arpack.h HINTS ARPACK_INCLUDE_DIRECTORY ${_quiet})
    if (ARPACK_LIBRARIES)
      add_library(ARPACK::ARPACK INTERFACE IMPORTED)
      target_link_libraries(ARPACK::ARPACK INTERFACE ${ARPACK_LIBRARIES})
      set(_arpack_found_method "library")
      set(_arpack_found_msg "as library ${ARPACK_LIBRARIES}")
      if (ARPACK_INCLUDE_DIRECTORIES)
        set_target_properties(
          ARPACK::ARPACK PROPERTIES
          INCLUDE_DIRECTORIES ${ARPACK_INCLUDE_DIRECTORIES}
        )
        set(
          _arpack_found_msg
          "${_arpack_found_msg} with include directory ${ARPACK_INCLUDE_DIRECTORIES}"
        )
      endif ()
    endif ()
  endif ()

  if (TARGET ARPACK::ARPACK)
    set(ARPACK_FOUND TRUE)

    # Find and link dependency if not found via CMake export file
    if (NOT _arpack_found_method STREQUAL "cmake")
      if (NOT TARGET BLAS::BLAS)
        find_package(BLAS REQUIRED ${_quiet})
      endif ()
      if (NOT TARGET LAPACK::LAPACK)
        find_package(LAPACK REQUIRED ${_quiet})
      endif ()
      target_link_libraries(ARPACK::ARPACK INTERFACE LAPACK::LAPACK BLAS::BLAS)
    endif ()
  endif ()

  unset(_quiet)
  unset(_required)

endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CustomArpack DEFAULT_MSG _arpack_found_msg)

unset(_arpack_found_method)
unset(_arpack_found_msg)
