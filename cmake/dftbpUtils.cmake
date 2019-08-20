# Replaces the extension of a given file
#
# Args:
#     oldext [in]: Old extension
#     newext [in]: New extension
#     fname [in]: File name in which extension should be replaced.
#     newfname [out]: File name after extension replacement.
#
function(dftbp_replace_extension oldext newext fname newfname)

  string(REGEX REPLACE "\\.${oldext}$" ".${newext}" _newfname ${fname})
  set(${newfname} ${_newfname} PARENT_SCOPE)

endfunction()


# Registers files for preprocessing
#
# Args:
#     preproc [in]: Preprocessor to use
#     preprocopts [in]:  Preprocessor command line arguments (but not in/out file)
#     oldext [in]: Extension of the unpreprocessed files.
#     newext [in]: Extension of the preprocessed files.
#     oldfiles [in]: List of unpreprocessed file names.
#     newfiles [out]: List of preprocessed file names.
#
function(dftbp_preprocess preproc preprocopts oldext newext oldfiles newfiles)

  set(_newfiles)
  foreach(oldfile IN LISTS oldfiles)
    dftbp_replace_extension(${oldext} ${newext} ${oldfile} newfile)
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${newfile}
      COMMAND ${preproc} ${preprocopts} ${CMAKE_CURRENT_SOURCE_DIR}/${oldfile} ${CMAKE_CURRENT_BINARY_DIR}/${newfile}
      MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${oldfile})
    list(APPEND _newfiles ${CMAKE_CURRENT_BINARY_DIR}/${newfile})
  endforeach()
  set(${newfiles} ${_newfiles} PARENT_SCOPE)

endfunction()


# Build -D command line arguments for Fypp preprocessor based on current configuration
#
# Args:
#     fyppflags [inout]: Current Fypp flags on enter, with -D options extended flags on exit.
#
function (dftbp_add_fypp_defines fyppflags)

  set(_fyppflags "${${fyppflags}}")
  if(WITH_OMP)
    list(APPEND _fyppflags -DWITH_OMP)
  endif()

  if(WITH_ARPACK)
    list(APPEND _fyppflags -DWITH_ARPACK)
  endif()

  if(WITH_DFTD3)
    list(APPEND _fyppflags -DWITH_DFTD3)
  endif()

  if(WITH_MPI)
    list(APPEND _fyppflags -DWITH_MPI -DWITH_SCALAPACK)
  endif()

  if(WITH_SOCKETS)
    list(APPEND _fyppflags -DWITH_SOCKETS)
  endif()

  if(WITH_ELSI)
    list(APPEND _fyppflags -DWITH_ELSI)
  endif()

  if(WITH_PEXSI)
    list(APPEND _fyppflags -DWITH_PEXSI)
  endif()

  if(WITH_GPU)
    list(APPEND _fyppflags -DWITH_GPU)
  endif()

  if(WITH_TRANSPORT)
    list(APPEND _fyppflags -DWITH_TRANSPORT)
  endif()

  set(${fyppflags} ${_fyppflags} PARENT_SCOPE)

endfunction()


# Gets DFTB+ release information.
#
# Args:
#   release [out]: Release string.
#
function(dftbp_get_release_name release)

  if(NOT EXISTS ${CMAKE_BINARY_DIR}/RELEASE)
    if(EXISTS ${CMAKE_SOURCE_DIR}/RELEASE)
      file(COPY ${CMAKE_SOURCE_DIR}/RELEASE DESTINATION ${CMAKE_BINARY_DIR})
    else()
      execute_process(
	COMMAND ${CMAKE_SOURCE_DIR}/utils/build/update_release ${CMAKE_BINARY_DIR}/RELEASE
	RESULT_VARIABLE exitcode)
      if(NOT exitcode EQUAL 0)
	file(WRITE ${CMAKE_BINARY_DIR}/RELEASE "(UNKNOWN RELEASE)")
      endif()
    endif()
  endif()
  file(READ ${CMAKE_BINARY_DIR}/RELEASE _release)
  separate_arguments(_release)
  set(${release} "${_release}" PARENT_SCOPE)

endfunction()


# Finds libraries and turns them into imported library targets
#
# Args:
#     libraries [in]: List of the library names to look for.
#     libpaths [in]: List of the paths to be used as hints when looking for the libraries.
#
function (dftbp_create_library_targets libraries libpaths)
  foreach(lib IN LISTS libraries)
    if(NOT (TARGET ${lib}))
      find_library(LIBPATH_FOUND ${lib} PATHS ${libpaths})
      IF(LIBPATH_FOUND)
	message(STATUS "Found library ${LIBPATH_FOUND}")
        add_library(${lib} IMPORTED UNKNOWN)
        set_target_properties(${lib} PROPERTIES IMPORTED_LOCATION ${LIBPATH_FOUND})
      else()
        message(FATAL_ERROR
	  "Could not find library '${lib}' using library path hints '${libpaths}'")
      endif()
      unset(LIBPATH_FOUND CACHE)
    endif()
  endforeach()
endfunction()


# Converts a space separated string into a list.
#
# Args:
#     * [in]: Name of the variables to convert. On exit the variables contain
#         the converted values.
#
function (dftbp_convert_to_list)
  foreach(varname IN LISTS ARGN)
    set(buffer "${${varname}}")
    separate_arguments(buffer)
    set(${varname} "${buffer}" PARENT_SCOPE)
  endforeach()
endfunction()


# Checks the build configuration on consistency and stops in case of inconsistencies
function (dftbp_ensure_config_consistency)

  if(WITH_ELSI AND NOT WITH_MPI)
    message(FATAL_ERROR "Buliding with ELSI requires MPI-parallel build enabled")
  endif()

  if(WITH_PEXSI AND (NOT WITH_MPI OR NOT WITH_ELSI))
    message(FATAL_ERROR "Building with PEXSI requires MPI-parallel build and ELSI enabled")
  endif()

  if(WITH_ARPACK AND WITH_MPI)
    message(FATAL_ERROR "Building with ARPACK requires MPI-parallel build disabled")
  endif()

  #if(BUILD_API AND WITH_MPI)
  #  message(FATAL_ERROR "Currently API can only be built if MPI-parallel build is disabled")
  #endif()

  if(("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG") AND WITH_OMP)
    message(FATAL_ERROR
      "NAG compiler usually creates crashing binary with OpenMP-parallelisation in debug mode. \
Disable OpenMP (WITH_OMP) when compiling in debug mode")
  endif()

  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    string(FIND "${CMAKE_Fortran_FLAGS}" "-standard-semantics" pos1)
    string(FIND "${CMAKE_Fortran_FLAGS}" "realloc_lhs" pos2)
    string(FIND "${CMAKE_Fortran_FLAGS}" "norealloc_lhs" pos3)
    if(NOT ((NOT pos1 EQUAL -1) OR ((NOT pos2 EQUAL -1) AND (pos3 EQUAL -1))))
      message(FATAL_ERROR "Intel compiler needs either the '-standard-semantics' or the '-assume \
realloc_lhs' option to produce correctly behaving (Fortran standard complying) code")
    endif()
  endif()

endfunction()


# Prepends a given path to a list of iles
#
# Args:
#     result [out]: Variable containing the results
#     prefix [in]: Prefix to add to each item
#     * [in]: List of items to be prefixed
#
function(dftbp_prepend_path result path)
  set(_result)
  foreach(var IN LISTS ARGN)
    list(APPEND _result "${path}/${var}")
  endforeach()
   set(${result} "${_result}" PARENT_SCOPE)
endfunction()
