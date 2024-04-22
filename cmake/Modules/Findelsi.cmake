# Distributed under the OSI-approved BSD 3-Clause License.
#
# Copyright (C) 2024  DFTB+ developers group
#

#[=======================================================================[.rst:
Findelsi
--------

Finds the ELSI library


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``elsi::elsi``
  The ELSI library

And, depending on the library build options:

``elsi::pexsi``
  The PEXSI library

Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variable:

``ELSI_FOUND``
  True if the system has the ELSI library

``ELSI_VERSION``
  Detected version of the library

#]=======================================================================]

if(NOT TARGET elsi::elsi)
  find_package(PkgConfig QUIET)
  pkg_check_modules(ELSI QUIET elsi)
  if(ELSI_FOUND)
    message(STATUS "Found 'elsi' via pkg-config")

    add_library(elsi::elsi INTERFACE IMPORTED)
    target_link_libraries(
      elsi::elsi
      INTERFACE
      "${ELSI_LINK_LIBRARIES}"
    )
    target_include_directories(
      elsi::elsi
      INTERFACE
      "${ELSI_INCLUDE_DIRS}"
    )
    target_include_directories(
      elsi::elsi
      INTERFACE
      "${ELSI_INCLUDEDIR}"
    )

    foreach(_lib IN LISTS ELSI_LINK_LIBRARIES)
      if(_lib MATCHES "pexsi")
        add_library(elsi::pexsi IMPORTED STATIC)
        set_target_properties(elsi::pexsi PROPERTIES
            IMPORTED_LOCATION ${_lib}
            IMPORTED_LINK_INTERFACE_LANGUAGES "CXX;Fortran")
        target_link_libraries(elsi::elsi INTERFACE elsi::pexsi)
        break()
      endif()
    endforeach()
    unset(_lib)

    # libmbd checks for lowercase variable name
    set(elsi_FOUND ${ELSI_FOUND})

    # DFTB+ checks for the lowercase variable name
    set(elsi_VERSION "${ELSI_VERSION}")
  endif()
endif()
