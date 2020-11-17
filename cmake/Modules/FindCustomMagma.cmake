# Distributed under the OSI-approved BSD 3-Clause License.
#
# Copyright (C) 2020 DFTB+ developers group
#

#[=======================================================================[.rst:
FindCustomMagma
---------------

Finds the MAGMA library


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``Magma::Magma``
  The MAGMA library


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variable:

``MAGMA_FOUND``
  True if the system has the MAGMA library


Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:

``MAGMA_DETECTION``
  Whether MAGMA libraries should be detected (default: True). If set to False,
  the settings in ``MAGMA_LIBRARY`` and ``MAGMA_INCLUDE_DIR`` will be used without any further
  checks.

``MAGMA_LIBRARY``

  Customized MAGMA library/libraries to use (instead of autodetected one). If
  not set or empty, the default library (magma) will be tried. The listed
  libraries will be checked for existence (unless disabled in
  ``MAGMA_DETECTION``) and the variable is overwritten to contain the libraries
  with their with full path.

``MAGMA_LIBRARY_DIR``
  Directories which should be looked up in order to find the customized libraries.

``MAGMA_INCLUDE_DIRECTORY``
  Customized MAGMA include directory.  

#]=======================================================================]

include(FindPackageHandleStandardArgs)
include(CustomLibraryFinder)

if(TARGET Magma::Magma)

  set(CUSTOMMAGMA_FOUND True)
  set(CustomMagma_FOUND True)
  set(MAGMA_FOUND True)
  set(Magma_FOUND True)

else()

  option(MAGMA_DETECTION "Whether MAGMA library should be detected" TRUE)

  if(MAGMA_DETECTION)

    find_package(PkgConfig QUIET)
    pkg_check_modules(_magma QUIET magma)
    # Overwrite PkgConfig values by user defined input if present.
    if(NOT "${MAGMA_LIBRARY}" STREQUAL "")
      set(_magma_LIBRARIES ${MAGMA_LIBRARY})
      set(_magma_LIBRARY_DIRS ${MAGMA_LIBRARY_DIR})
    endif()
    if(NOT "${MAGMA_INCLUDE_DIR}" STREQUAL "")
      set(_magma_INCLUDE_DIRS ${MAGMA_INCLUDE_DIR})
    endif()

    find_custom_libraries("${_magma_LIBRARIES}" "${_magma_LIBRARY_DIRS}"
      "${CustomMagma_FIND_QUIETLY}" _libs)
    set(MAGMA_LIBRARY "${_libs}" CACHE STRING "List of MAGMA libraries to link" FORCE)
    unset(_libs)
    unset(_magma_LIBRARIES)
    unset(_magma_LIBRARY_DIRS)

    # Check include file
    find_path(MAGMA_INCLUDE_DIRECTORY NAMES magma.mod PATHS ${_magma_INCLUDE_DIRS}
      PATH_SUFFIXES magma)
    unset(_magma_INCLUDE_DIRS)
  
    set(MAGMA_DETECTION FALSE CACHE BOOL "Whether MAGMA library should be detected" FORCE)
    
  endif()

  find_package_handle_standard_args(CustomMagma REQUIRED_VARS MAGMA_LIBRARY MAGMA_INCLUDE_DIRECTORY)

  set(MAGMA_FOUND ${CUSTOMMAGMA_FOUND})
  set(Magma_FOUND ${CUSTOMMAGMA_FOUND})

  if(MAGMA_FOUND AND NOT TARGET Magma::Magma)
    add_library(Magma::Magma INTERFACE IMPORTED)
    target_link_libraries(Magma::Magma INTERFACE "${MAGMA_LIBRARY}")
    target_include_directories(Magma::Magma INTERFACE "${MAGMA_INCLUDE_DIRECTORY}")
  endif()

endif()
