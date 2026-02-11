# Distributed under the OSI-approved BSD 3-Clause License.
#
# Copyright (C) 2025  DFTB+ developers group
#

#[=======================================================================[.rst:
Findelpa
--------

Finds the ELPA library


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``elpa::elpa``
  The ELPA library

Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variable:

``ELPA_FOUND``
  True if the system has the ELPA library

``ELPA_VERSION``
  Detected version of the library

#]=======================================================================]

if(NOT TARGET elpa::elpa)
  find_package(PkgConfig QUIET)
  pkg_check_modules(ELPA QUIET elpa)
  if(ELPA_FOUND)
    message(STATUS "Found 'elpa' via pkg-config")

    add_library(elpa::elpa INTERFACE IMPORTED)
    target_link_libraries(
      elpa::elpa
      INTERFACE
      "${ELPA_LINK_LIBRARIES}"
    )
    target_include_directories(
      elpa::elpa
      INTERFACE
      "${ELPA_INCLUDE_DIRS}/modules"
    )

    # DFTB+ checks for the lowercase variable name
    set(elpa_VERSION "${ELPA_VERSION}")
  endif()
endif()
