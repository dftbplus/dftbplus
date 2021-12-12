# Distributed under the OSI-approved BSD 2-Clause License.
#
# Copyright (C) 2020 DFTB+ developers group
#

#[=======================================================================[.rst:
FindCustomScalapack
-------------------

Finds the ScaLAPACK library.

This is a simple auto-detection module for the ScaLAPACK library (looks
basically for a library with the name 'scalapack'). It also assumes that
ScaLAPACK is MPI-based and defines a respective dependency on
``MPI::MPI_Fortran.``


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``scalapack``
  The ScaLAPACK library. The name was chosen to be compatible with the name used by the CMake
  export file (which is unfortunately broken in several binary packages).

``Scalapack::Scalapack``
  Name space (and CMake-style capitalized) variant of the scalapack target.


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variable:

``SCALAPACK_FOUND``
  True if the system has the SCALAPACK library


Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``SCALAPACK_DETECTION``
  Whether ScaLAPACK libraries should be detected (default: True). If set to False,
  the settings in ``SCALAPACK_LIBRARY`` will be used without any further checks.

``SCALAPACK_LIBRARY``
  Customized ScaLAPACK library/libraries to use.  If no SCALAPACK library is
  required (e.g. the linker automatically links it) set
  ``SCALAPACK_LIBRARY="NONE"``. If not set or empty, it will use
  ``find_package(scalapack)`` to find the scalapack library (relying on the
  existence of a CMake export file for scalapack). Otherwise, the listed
  libraries will be checked for existence (unless disabled in
  ``SCALAPACK_DETECTION``) and the variable is overwritten to contain the
  libraries with their with full path.

  Note: On Ubuntu Focal (20.4 LTS) system, the installed CMake export file is broken.
  You should not realy on the CMake-autodetection on those systems.

``SCALAPACK_LIBRARY_DIR``
  Directories which should be looked up in order to find the customized libraries.

#]=======================================================================]

include(FindPackageHandleStandardArgs)
include(CustomLibraryFinder)

if(TARGET scalapack)

  set(CUSTOMSCALAPACK_FOUND True)
  set(CustomScalapack_FOUND True)
  set(SCALAPACK_FOUND True)
  set(Scalapack_FOUND True)

else()

  find_package(MPI ${Find_Scalapack_REQUIRED})

  option(SCALAPACK_DETECTION "Whether ScaLAPACK library should be detected" TRUE)

  if(SCALAPACK_DETECTION)

    if("${SCALAPACK_LIBRARY}" STREQUAL "")

      # Try Scalapack via CMake export file
      find_package(scalapack QUIET)
      if(scalapack_FOUND)
        get_target_property(_scalapack_library scalapack IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG)
        set(SCALAPACK_LIBRARY "${_scalapack_library}" CACHE STRING "ScaLAPACK library to link"
          FORCE)
        unset(_scalapack_library)
      else()
        # Very simple ScaLAPACK auto-detection: looking for a library called scalapack
        find_library(SCALAPACK_LIBRARY scalapack HINTS ${SCALAPACK_LIBRARY_DIR})
      endif()

    elseif(NOT "${SCALAPACK_LIBRARY}" STREQUAL "NONE")

      find_custom_libraries("${SCALAPACK_LIBRARY}" "${SCALAPACK_LIBRARY_DIR}"
        "${CustomScalapack_FIND_QUIETLY}" _libs)
      set(SCALAPACK_LIBRARY "${_libs}" CACHE STRING "List of ScaLAPACK libraries to link" FORCE)
      unset(_libs)

    endif()

    set(SCALAPACK_DETECTION False CACHE BOOL "Whether ScaLAPACK libraries should be detected" FORCE)

  endif()

  find_package_handle_standard_args(CustomScalapack REQUIRED_VARS SCALAPACK_LIBRARY
    MPI_Fortran_FOUND)

  set(SCALAPACK_FOUND ${CUSTOMSCALAPACK_FOUND})
  set(Scalapack_FOUND ${CUSTOMSCALAPACK_FOUND})

  if(SCALAPACK_FOUND)
    if(NOT TARGET scalapack)
      add_library(scalapack INTERFACE IMPORTED)
      if(NOT "${SCALAPACK_LIBRARY}" STREQUAL "NONE")
        target_link_libraries(scalapack INTERFACE "${SCALAPACK_LIBRARY}")
      endif()
      target_link_Libraries(scalapack INTERFACE MPI::MPI_Fortran)
    endif()
  endif()

  mark_as_advanced(SCALAPACK_DETECTION SCALAPACK_LIBRARY SCALAPACK_LIBRARY_DIR)

endif()

# Add namespaced library name variant
if(TARGET scalapack AND NOT TARGET Scalapack::Scalapack)
  add_library(Scalapack::Scalapack INTERFACE IMPORTED)
  target_link_libraries(Scalapack::Scalapack INTERFACE scalapack)
endif()
