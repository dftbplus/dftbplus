# Simple minded finder for MAGMA
#
# Defines following variables:
#
#   MAGMA_FOUND
#   MAGMA_LIBRARIES
#   MAGMA_LIBRARY_DIRS
#   MAGMA_INCLUDE_DIRS
#
# provided the variables had not been defined already. If that is the case, they
# will be passed unaltered.
#
if(NOT DEFINED MAGMA_FOUND)
  find_package(PkgConfig)

  if (NOT DEFINED MAGMA_LIBRARIES OR MAGMA_LIBRARIES STREQUAL "")
    pkg_check_modules(_pc_magma QUIET magma)
    set(_magma_LIBRARIES ${_pc_magma_LIBRARIES})
    set(_magma_LIBRARY_DIRS ${_pc_magma_LIBRARY_DIRS})
  else()
    set(_magma_LIBRARIES ${MAGMA_LIBRARIES})
    set(_magma_LIBRARY_DIRS ${MAGMA_LIBRARY_DIRS})
  endif()

  if (NOT DEFINED MAGMA_INCLUDE_DIRS OR MAGMA_INCLUDE_DIRS STREQUAL "")
    find_path(_magma_INCLUDE_DIR NAMES magma.mod QUIET)
  else()
    # Consistency check: whether provided path contains "magma.mod"
    find_file(_magma_mod NAMES magma.mod PATHS ${MAGMA_INCLUDE_DIRS})
    if(_magma_mod)
      get_filename_component(_magma_INCLUDE_DIR ${_magma_mod} DIRECTORY)
    endif()
  endif()

  mark_as_advanced(MAGMA_FOUND _magma_INCLUDE_DIR _magma_LIBRARIES _magma_LIBRARY_DIRS)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MAGMA
    FOUND_VAR MAGMA_FOUND
    REQUIRED_VARS _magma_LIBRARIES _magma_INCLUDE_DIR)
  
  if(MAGMA_FOUND)
    if(NOT DEFINED MAGMA_LIBRARIES OR MAGMA_LIBRARIES STREQUAL "")
      set(MAGMA_LIBRARIES ${_magma_LIBRARIES})
      set(MAGMA_LIBRARY_DIRS ${_magma_LIBRARY_DIRS})
    endif()
    if(NOT DEFINED MAGMA_INCLUDE_DIRS OR MAGMA_INCLUDE_DIRS STREQUAL "")
      set(MAGMA_INCLUDE_DIRS ${_magma_INCLUDE_DIR})
    endif()
  endif()

endif()
