# Distributed under the OSI-approved BSD 3-Clause License.
#
# Copyright (C) 2020 DFTB+ developers group
#

#[=======================================================================[.rst:
FindCustomPlumed
----------------

Finds the PLUMED library


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``Plumed::Plumed``
  The PLUMED library


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variable:

``PLUMED_FOUND``
  True if the system has the PLUMED library


Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``PLUMED_DETECTION``
  Whether PLUMED libraries should be detected (default: True). If set to False,
  the settings in ``PLUMED_LIBRARY`` will be used without any further
  checks.

``PLUMED_LIBRARY``

  Customized PLUMED library/libraries to use (instead of autodetected one). If
  not set or empty, the default library (plumed) will be tried. The listed
  libraries will be checked for existence (unless disabled in
  ``PLUMED_DETECTION``) and the variable is overwritten to contain the libraries
  with their with full path.

``PLUMED_LIBRARY_DIR``
  Directories which should be looked up in order to find the customized libraries.


#]=======================================================================]

include(FindPackageHandleStandardArgs)
include(CustomLibraryFinder)

if(TARGET Plumed::Plumed)

  set(CUSTOMPLUMED_FOUND True)
  set(CustomPlumed_FOUND True)
  set(PLUMED_FOUND True)
  set(Plumed_FOUND True)

else()

  option(PLUMED_DETECTION "Whether PLUMED library should be detected" TRUE)

  if(PLUMED_DETECTION)

    find_package(PkgConfig)
    pkg_check_modules(_plumed QUIET plumed)

    # Overwrite PkgConfig values by user defined input if present.
    if(NOT "${PLUMED_LIBRARY}" STREQUAL "")
      set(_plumed_LIBRARIES ${PLUMED_LIBRARY})
      set(_plumed_LIBRARY_DIRS ${PLUMED_LIBRARY_DIR})
    endif()

    find_custom_libraries("${_plumed_LIBRARIES}" "${_plumed_LIBRARY_DIRS}"
      "${CustomPlumed_FIND_QUIETLY}" _libs)
    set(PLUMED_LIBRARY "${_libs}" CACHE STRING "List of PLUMED libraries to link" FORCE)
    unset(_libs)
    unset(_plumed_LIBRARIES)
    unset(_plumed_LIBRARY_DIRS)

    set(PLUMED_DETECTION False CACHE BOOL "Whether PLUMED libraries should be detected" FORCE)

  endif()

  find_package_handle_standard_args(CustomPlumed REQUIRED_VARS PLUMED_LIBRARY)

  set(PLUMED_FOUND ${CUSTOMPLUMED_FOUND})
  set(Plumed_FOUND ${CUSTOMPLUMED_FOUND})

  if(PLUMED_FOUND AND NOT TARGET Plumed::Plumed)
    add_library(Plumed::Plumed INTERFACE IMPORTED)
    target_link_libraries(Plumed::Plumed INTERFACE "${PLUMED_LIBRARY}")
  endif()

endif()
