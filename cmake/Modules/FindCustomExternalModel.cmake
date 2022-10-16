# Distributed under the OSI-approved BSD 3-Clause License.
#
# Copyright (C) 2022  DFTB+ developers group
#

#[=======================================================================[.rst:
FindCustomExternalModel
-----------------------

Finds library with an external model


Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``ExternalModel::ExternalModel``
  The external model library


Result Variables
^^^^^^^^^^^^^^^^

This module will define the following variable:

``EXTERNALMODEL_FOUND``
  True if the system has the EXTERNALMODEL library


Cache variables
^^^^^^^^^^^^^^^

The following cache variables may be set to influence the library detection:


``EXTERNALMODEL_LIBRARY``

  Customized EXTERNALMODEL library/libraries to use (instead of
  autodetected one). If not set or empty, the default library
  (externalmodel) will be tried. The listed libraries will be checked
  for existence and the variable is overwritten to contain the
  libraries with their with full path.

``EXTERNALMODEL_LIBRARY_DIR``
  Directories which should be looked up in order to find the customized libraries.


#]=======================================================================]

include(FindPackageHandleStandardArgs)
include(CustomLibraryFinder)

if(TARGET ExternalModel::ExternalModel)

  set(CUSTOMEXTERNALMODEL_FOUND True)
  set(CustomExternalModel_FOUND True)
  set(EXTERNALMODEL_FOUND True)
  set(ExternalModel_FOUND True)

else()

  option(EXTERNALMODEL_DETECTION "Whether EXTERNALMODEL library should be detected" TRUE)

  if(EXTERNALMODEL_DETECTION)

    find_package(PkgConfig)
    pkg_check_modules(_externalmodel QUIET externalmodel)

    # Overwrite PkgConfig values by user defined input if present.
    if(NOT "${EXTERNALMODEL_LIBRARY}" STREQUAL "")
      set(_externalmodel_LIBRARIES ${EXTERNALMODEL_LIBRARY})
      set(_externalmodel_LIBRARY_DIRS ${EXTERNALMODEL_LIBRARY_DIR})
    endif()

    find_custom_libraries("${_externalmodel_LIBRARIES}" "${_externalmodel_LIBRARY_DIRS}"
      "${CustomExternalModel_FIND_QUIETLY}" _libs)
    set(EXTERNALMODEL_LIBRARY "${_libs}" CACHE STRING "List of EXTERNALMODEL libraries to link" FORCE)
    unset(_libs)
    unset(_externalmodel_LIBRARIES)
    unset(_externalmodel_LIBRARY_DIRS)

    set(EXTERNALMODEL_DETECTION False CACHE BOOL "Whether EXTERNALMODEL libraries should be detected" FORCE)

  endif()

  find_package_handle_standard_args(CustomExternalModel REQUIRED_VARS EXTERNALMODEL_LIBRARY)

  set(EXTERNALMODEL_FOUND ${CUSTOMEXTERNALMODEL_FOUND})
  set(ExternalModel_FOUND ${CUSTOMEXTERNALMODEL_FOUND})

  if(EXTERNALMODEL_FOUND AND NOT TARGET ExternalModel::ExternalModel)
    add_library(ExternalModel::ExternalModel INTERFACE IMPORTED)
    target_link_libraries(ExternalModel::ExternalModel INTERFACE "${EXTERNALMODEL_LIBRARY}")
  endif()

endif()
