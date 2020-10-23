# Distributed under the OSI-approved BSD 3-Clause License.
#
# Copyright (C) 2020 DFTB+ developers group
#

#[=======================================================================[.rst:
FindCustomArpack
----------------

Finds the ARPACK library


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``Arpack::Arpack``
  The ARPACK library


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variable:

``ARPACK_FOUND``
  True if the system has the ARPACK library


Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``ARPACK_DETECTION``
  Whether ARPACK libraries should be detected (default: True). If set to False,
  the settings in ``ARPACK_LIBRARY`` will be used without any further checks.

``ARPACK_LIBRARY``

  Customized ARPACK library/libraries to use (instead of autodetected one). If
  not set or empty, the default library (arpack) will be tried. The listed
  libraries will be checked for existence (unless disabled in
  ``ARPACK_DETECTION``) and the variable is overwritten to contain the libraries
  with their with full path.

``ARPACK_LIBRARY_DIR``
  Directories which should be looked up in order to find the customized libraries.

#]=======================================================================]

include(FindPackageHandleStandardArgs)
include(CustomLibraryFinder)

if(TARGET Arpack::Arpack)

  set(CUSTOMARPACK_FOUND True)
  set(CustomArpack_FOUND True)
  set(ARPACK_FOUND True)
  set(Arpack_FOUND True)

else()

  option(ARPACK_DETECTION "Whether ARPACK library should be detected" TRUE)

  if(ARPACK_DETECTION)

    if("${ARPACK_LIBRARY}" STREQUAL "")
      # Very simple ARPACK auto-detection
      find_library(ARPACK_LIBRARY arpack HINTS ${ARPACK_LIBRARY_DIR})

    else()

      # Library explicitely set by the user, search for those libraries
      find_custom_libraries("${ARPACK_LIBRARY}" "${ARPACK_LIBRARY_DIR}"
        "${CustomArpack_FIND_QUIETLY}" _libs)
      set(ARPACK_LIBRARY "${_libs}" CACHE STRING "List of ARPACK libraries to link" FORCE)
      unset(_libs)

    endif()

    set(ARPACK_DETECTION False CACHE BOOL "Whether ARPACK libraries should be detected" FORCE)

  endif()

  find_package_handle_standard_args(CustomArpack REQUIRED_VARS ARPACK_LIBRARY)

  set(ARPACK_FOUND ${CUSTOMARPACK_FOUND})
  set(Arpack_FOUND ${CUSTOMARPACK_FOUND})

  if(ARPACK_FOUND AND NOT TARGET Arpack::Arpack)
    add_library(Arpack::Arpack INTERFACE IMPORTED)
    target_link_libraries(Arpack::Arpack INTERFACE "${ARPACK_LIBRARY}")
  endif()

  mark_as_advanced(ARPACK_DETECTION ARPACK_LIBRARY ARPACK_LIBRARY_DIR)

endif()
