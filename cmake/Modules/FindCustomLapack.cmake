# Distributed under the OSI-approved BSD 2-Clause License.
#
# Copyright (C) 2020 DFTB+ developers group
#

#[=======================================================================[.rst:
FindCustomLapack
----------------

Finds the LAPACK library

This is a wrapper around CMakes FindLAPACK module with the additional
possibility to customize the library name manually. In latter case the module will
check the existence of those libraries and stop if they are not found.

Note: The module is named FindLapack (and not FindLAPACK) to avoid name
collision with CMakes built-in FindLAPACK module.


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``LAPACK::LAPACK``
  The LAPACK library


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variable:

``LAPACK_FOUND``
  True if the system has the LAPACK library


Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``LAPACK_DETECTION``
  Whether LAPACK libraries should be detected (default: True). If set to False,
  the settings in ``LAPACK_LIBRARY`` will be used without any further checks.

``LAPACK_LIBRARY``
  Customized LAPACK library/libraries to use (instead of autodetected ones).  If
  no LAPACK library is required (e.g. the linker automatically links it) set
  ``LAPACK_LIBRARY="NONE"``. If not set or empty, the built-in LAPACK finder
  (the findLAPACK module) will be invoked. Otherwise, the listed libraries
  will be checked for existence (unless disabled in ``LAPACK_DETECTION``) and
  the variable is overwritten to contain the libraries with their with full
  path.

``LAPACK_LIBRARY_DIR``
  Directories which should be looked up in order to find the customized libraries.

``LAPACK_LINKER_FLAG``
  Flags to use when linking LAPACK

Additionally, the cache variables of the built-in FindLAPACK modules may used to
influence the LAPACK detection if the built-in module is invoked.

#]=======================================================================]

include(FindPackageHandleStandardArgs)
include(CustomLibraryFinder)

if(TARGET LAPACK::LAPACK)

  set(CUSTOMLAPACK_FOUND True)
  set(CustomLapack_FOUND True)
  set(LAPACK_FOUND True)
  set(Lapack_FOUND True)

else()

  option(LAPACK_DETECTION "Whether LAPACK library should be detected" TRUE)

  if(LAPACK_DETECTION)
    # LAPACK has either not been found yet or it was found by an older built-in findLAPACK module.
    # which does not provide the imported target LAPACK::LAPACK

    if("${LAPACK_LIBRARY}" STREQUAL "")

      # No user customized LAPACK library, try built-in finder
      if(NOT LAPACK_FOUND)
        find_package(LAPACK)
      endif()
      set(LAPACK_LIBRARY "${LAPACK_LIBRARIES}" CACHE STRING "LAPACK library to link" FORCE)
      set(LAPACK_LINKER_FLAG "${LAPACK_LINKER_FLAGS}" CACHE STRING
        "Linker flags to use when linking LAPACK" FORCE)

    elseif(NOT "${LAPACK_LIBRARY}" STREQUAL "NONE")

      # LAPACK explicitely set by the user, search for those libraries
      find_custom_libraries("${LAPACK_LIBRARY}" "${LAPACK_LIBRARY_DIR}"
        "${CustomLapack_FIND_QUIETLY}" _libs)
      set(LAPACK_LIBRARY "${_libs}" CACHE STRING "List of LAPACK libraries to link" FORCE)
      unset(_libs)

    endif()

    set(LAPACK_DETECTION False CACHE BOOL "Whether LAPACK libraries should be detected" FORCE)

  endif()

  find_package_handle_standard_args(CustomLapack REQUIRED_VARS LAPACK_LIBRARY)

  set(LAPACK_FOUND ${CUSTOMLAPACK_FOUND})
  set(Lapack_FOUND ${CUSTOMLAPACK_FOUND})


  if(LAPACK_FOUND AND NOT TARGET LAPACK::LAPACK)

    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    if(NOT "${LAPACK_LIBRARY}" STREQUAL "NONE")
      target_link_libraries(LAPACK::LAPACK INTERFACE "${LAPACK_LIBRARY}")
    endif()
    if(NOT "${LAPACK_LINKER_FLAG}" STREQUAL "")
      target_link_options(LAPACK::LAPACK INTERFACE "${LAPACK_LINKER_FLAG}")
    endif()
  endif()

  mark_as_advanced(LAPACK_DETECTION LAPACK_LIBRARY LAPACK_LIBRARY_DIR LAPACK_LINKER_FLAG)

endif()
