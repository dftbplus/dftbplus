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

  string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
  if("${CMAKE_BUILD_TYPE_UPPER}" STREQUAL "DEBUG")
    list(APPEND _fyppflags -DDEBUG=1)
  else()
    list(APPEND _fyppflags -DDEBUG=0)
  endif()

  if(INTERNAL_ERFC)
    list(APPEND _fyppflags -DINTERNAL_ERFC=1)
  endif()

  if(WITH_OMP)
    list(APPEND _fyppflags -DWITH_OMP)
  endif()

  if(WITH_ARPACK)
    list(APPEND _fyppflags -DWITH_ARPACK)
  endif()

  if(WITH_DFTD3)
    list(APPEND _fyppflags -DWITH_DFTD3)
  endif()

  if(WITH_PLUMED)
    list(APPEND _fyppflags -DWITH_PLUMED)
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

  if(ELSI_WITH_PEXSI)
    list(APPEND _fyppflags -DWITH_PEXSI)
  endif()

  if(WITH_GPU)
    list(APPEND _fyppflags -DWITH_GPU)
  endif()

  if(WITH_TRANSPORT)
    list(APPEND _fyppflags -DWITH_TRANSPORT)
  endif()

  if(WITH_C_EXECUTABLES)
    list(APPEND _fyppflags -DWITH_C_EXECUTABLES)
  endif()

  if(BUILD_SHARED_LIBS)
    list(APPEND _fyppflags -DBUILD_SHARED_LIBS)
  endif()

  set(${fyppflags} ${_fyppflags} PARENT_SCOPE)

endfunction()


# Gets DFTB+ release information.
#
# Args:
#   release [out]: Release string.
#
function(dftbp_get_release_name release)

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
  file(READ ${CMAKE_BINARY_DIR}/RELEASE _release)
  separate_arguments(_release)
  set(${release} "${_release}" PARENT_SCOPE)

endfunction()


# Gets DFTB+ API version information.
#
# Args:
#   apiversion [out]: Version string.
#   apimajor [out]: Major release number (as string).
#   apiminor [out]: Minor release number (as string).
#   apipatch [out]: Patch release number (as string).
#
function(dftbp_get_api_version apiversion apimajor apiminor apipatch)

  file(STRINGS ${CMAKE_SOURCE_DIR}/prog/dftb+/api/mm/API_VERSION _api
    REGEX "^[0-9]+\.[0-9]+\.[0-9]+$")
  string(REGEX MATCHALL "[0-9]+" _api_list "${_api}")
  list(GET _api_list 0 _api_major)
  list(GET _api_list 1 _api_minor)
  list(GET _api_list 2 _api_patch)

  set(${apiversion} "${_api}" PARENT_SCOPE)
  set(${apimajor} "${_api_major}" PARENT_SCOPE)
  set(${apiminor} "${_api_minor}" PARENT_SCOPE)
  set(${apipatch} "${_api_patch}" PARENT_SCOPE)

endfunction()


# Finds libraries and turns them into imported library targets
#
# Args:
#     libraries [in]: List of the library names to look for.
#     libpaths [in]: List of the paths to be used as hints when looking for the libraries.
#
function (dftbp_create_library_targets libraries libpaths)
  foreach(lib IN LISTS libraries)
    string(REGEX MATCH "^[ ]*-.*" option ${lib})
    # If the library is a linker option, skip target conversion (use it literally)
    if(NOT "${option}" STREQUAL "")
      continue()
    endif()
    if(TARGET ${lib})
      continue()
    endif()
    if(IS_ABSOLUTE "${lib}")
      if(EXISTS "${lib}")
        message(STATUS "Found library file ${lib}")
      else()
        message(FATAL_ERROR "Library file '${lib}' not found")
      endif()
    else()
      find_library(LIBPATH ${lib} HINTS ${libpaths})
      if(LIBPATH)
	message(STATUS "Found library ${lib} as ${LIBPATH}")
	if(EXPORT_EXTLIBS_WITH_PATH)
	  add_library(${lib} INTERFACE)
	  target_link_libraries(${lib} INTERFACE ${LIBPATH})
	  install(TARGETS ${lib} EXPORT dftbplus-targets)
	else()
          add_library(${lib} IMPORTED UNKNOWN)
          set_target_properties(${lib} PROPERTIES IMPORTED_LOCATION ${LIBPATH})
	endif()
      else()
        message(FATAL_ERROR
	  "Could not find library '${lib}' using library path hints '${libpaths}'")
      endif()
      unset(LIBPATH CACHE)
    endif()
  endforeach()
endfunction()


# Checks the build configuration on consistency and stops in case of inconsistencies
function (dftbp_ensure_config_consistency)

  if(WITH_ELSI AND NOT WITH_MPI)
    message(FATAL_ERROR "Building with ELSI requires MPI-parallel build enabled")
  endif()

  if(WITH_PEXSI AND (NOT WITH_MPI OR NOT WITH_ELSI))
    message(FATAL_ERROR "Building with PEXSI requires MPI-parallel build and ELSI enabled")
  endif()

  if(WITH_ARPACK AND WITH_MPI)
    message(FATAL_ERROR "Building with ARPACK requires MPI-parallel build disabled")
  endif()

  if(WITH_GPU AND WITH_MPI)
    message(FATAL_ERROR "Building with GPU support and MPI parallelisation disabled")
  endif()

  string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
  if(("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
      AND ("${CMAKE_BUILD_TYPE_UPPER}" STREQUAL "DEBUG") AND WITH_OMP)
    message(FATAL_ERROR
      "NAG compiler usually creates crashing binary with OpenMP-parallelisation in debug mode. \
Disable OpenMP (WITH_OMP) when compiling in debug mode")
  endif()

  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    string(FIND "${CMAKE_Fortran_FLAGS}" "-standard-semantics" pos1)
    string(FIND "${CMAKE_Fortran_FLAGS}" "realloc_lhs" pos2)
    string(FIND "${CMAKE_Fortran_FLAGS}" "norealloc_lhs" pos3)
    if(NOT ((NOT pos1 EQUAL -1) OR ((NOT pos2 EQUAL -1) AND (pos3 EQUAL -1))))
      message(FATAL_ERROR "Intel Fortran compiler needs either the '-standard-semantics' or the "
        "'-assume realloc_lhs' option to produce correctly behaving (Fortran standard complying) "
        "code")
    endif()
  endif()

  set(pkgconfig_languages C Fortran)
  list(FIND pkgconfig_languages "${PKGCONFIG_LANGUAGE}" pos)
  if(pos EQUAL -1)
    string(REPLACE ";" "\", \"" pkgconfig_languages_str "${pkgconfig_languages}")
    set(pkgconfig_languages_str "\"${pkgconfig_languages_str}\"")
    message(FATAL_ERROR
      "Invalid language \"${PKGCONFIG_LANGUAGE}\" for PKGCONFIG_LANGUAGE (possible choices: ${pkgconfig_languages_str})")
  endif()

endfunction()


# Prepends a given prefix to a list of items.
#
# Args:
#     result [out]: Variable containing the results
#     prefix [in]: Prefix to add to each item
#     * [in]: List of items to be prefixed
#
function(dftbp_add_prefix prefix itemlist result)
  set(_result)
  foreach(var IN LISTS itemlist)
    list(APPEND _result "${prefix}${var}")
  endforeach()
  set(${result} "${_result}" PARENT_SCOPE)
endfunction()


# Returns the parameters needed to create a pkg-config export file
#
# Args:
#     pkgconfig_requires [out]: Value for the Requires field.
#     pkgconfig_libs [out]: Value for the Libs field.
#     pkgconfig_libs_private [out]: Value for the Libs.private field.
#     pkgconfig_c_flags [out]: value for the cflags field.
#
function(dftbp_get_pkgconfig_params pkgconfig_requires pkgconfig_libs pkgconfig_libs_private
    pkgconfig_c_flags)

  set(_pkgconfig_libs "-L${INSTALL_LIB_DIR} -ldftbplus")
  set(_pkgconfig_libs_private)

  dftbp_add_prefix("-l" "${EXPORTED_COMPILED_LIBRARIES}" complibs)
  list(APPEND _pkgconfig_libs "${complibs}")

  dftbp_add_prefix("-L" "${EXPORTED_EXTERNAL_LIBRARY_DIRS}" extlibdirs)
  list(APPEND _pkgconfig_libs_private "${extlibdirs}")

  dftbp_library_linking_flags("${EXPORTED_EXTERNAL_LIBRARIES}" extlibs)
  list(APPEND _pkgconfig_libs_private "${extlibs}")

  if(PKGCONFIG_LANGUAGE STREQUAL "C")

    set(implibdirs "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
    list(REMOVE_ITEM implibdirs ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES})
    dftbp_add_prefix("-L" "${implibdirs}" implibdirs)
    list(APPEND _pkgconfig_libs_private "${implibdirs}")

    set(implibs "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
    list(REMOVE_ITEM implibs ${CMAKE_C_IMPLICIT_LINK_LIBRARIES})
    dftbp_library_linking_flags("${implibs}" implibs)
    list(APPEND _pkgconfig_libs_private "${implibs}")

    set(_pkgconfig_c_flags "-I${INSTALL_INC_DIR} ${CMAKE_C_FLAGS}")

  else()

    set(implibdirs "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
    list(REMOVE_ITEM implibdirs ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
    dftbp_add_prefix("-L" "${implibdirs}" implibdirs)
    list(APPEND _pkgconfig_libs_private "${implibdirs}")

    set(implibs "${CMAKE_C_IMPLICIT_LINK_LIBRARIES}")
    list(REMOVE_ITEM implibs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
    dftbp_library_linking_flags("${implibs}" implibs)
    list(APPEND _pkgconfig_libs_private "${implibs}")

    set(_pkgconfig_c_flags "-I${INSTALL_MOD_DIR} ${CMAKE_Fortran_FLAGS}")

  endif()

  string(REPLACE ";" " " _pkgconfig_libs "${_pkgconfig_libs}")
  string(REPLACE ";" " " _pkgconfig_libs_private "${_pkgconfig_libs_private}")

  set(_pkgconfig_libs_private "${_pkgconfig_libs_private} ${CMAKE_EXE_LINKER_FLAGS}")
  string(REPLACE ";" " " _pkgconfig_requires "${EXPORTED_EXTERNAL_PACKAGES}")

  set(${pkgconfig_requires} "${_pkgconfig_requires}" PARENT_SCOPE)
  set(${pkgconfig_libs} "${_pkgconfig_libs}" PARENT_SCOPE)
  set(${pkgconfig_libs_private} "${_pkgconfig_libs_private}" PARENT_SCOPE)
  set(${pkgconfig_c_flags} "${_pkgconfig_c_flags}" PARENT_SCOPE)

endfunction()


# Returns library linking flags for given libraries.
#
# If the library is a full path, it is linked directly, otherwise with the "-l" option
#
# Args:
#     libraries [in]: List of libraries to check.
#     linkflags [out]: List of flags to link the libraries
#
function(dftbp_library_linking_flags libraries linkflags)
  set(_linkflags)
  foreach(library IN LISTS libraries)
    string(FIND "${library}" "-" dashpos)
    if(dashpos EQUAL 1)
      list(APPEND _linkflags "${library}")
    elseif(IS_ABSOLUTE "${library}")
      list(APPEND _linkflags "${library}")
    else()
      list(APPEND _linkflags "-l${library}")
    endif()
  endforeach()
  set(${linkflags} "${_linkflags}" PARENT_SCOPE)
endfunction()


# Stops the code if the source and the build folders are identical.
#
function(dftbp_ensure_out_of_source_build)

  get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

  if("${srcdir}" STREQUAL "${bindir}")
    message(FATAL_ERROR
      "It is not allowed to configure and build DFTB+ from its source folder. Please, create a \
separate build directory and invoke CMake from that directory. See the INSTALL.rst file for \
detailed build instructions.")
  endif()

endfunction()


# Makes sure, that a compiler has been already defined for a given language
#
# Args:
#     language [in]: The language to look at.
#
function(dftbp_ensure_compiler_def language)

  if(NOT DEFINED CMAKE_${language}_COMPILER)
    message(FATAL_ERROR "Undefined ${language} compiler. The automatic detection of compilers, \
flags and libraries is disabled. You must provide configuration parameters explicitely (e.g. in a \
toolchain file). See the INSTALL.rst file for detailed instructions.")
  endif()

endfunction()


# Loads global build settings (either from config.cmake or from user defined file)
#
macro (dftbp_load_build_settings)

  if(NOT DEFINED BUILD_CONFIG_FILE)
    if(DEFINED ENV{DFTBPLUS_BUILD_CONFIG_FILE} AND NOT ENV{DFTBPLUS_BUILD_CONFIG_FILE} STREQUAL "")
      set(BUILD_CONFIG_FILE "$ENV{DFTBPLUS_BUILD_CONFIG_FILE}")
    else()
      set(BUILD_CONFIG_FILE "${CMAKE_SOURCE_DIR}/config.cmake")
    endif()
  endif()
  message(STATUS "Reading global build config file: ${BUILD_CONFIG_FILE}")
  include(${BUILD_CONFIG_FILE})
  
endmacro()


# Tries to guess which toolchain to load based on the environment.
#
# Args:
#     toolchain [out]: Name of the selected toolchain or undefined if it could not be selected
#
function(dftbp_guess_toolchain toolchain)

  if("${CMAKE_Fortran_COMPILER_ID}|${CMAKE_C_COMPILER_ID}" STREQUAL "GNU|GNU")
    set(_toolchain "gnu")
  elseif("${CMAKE_Fortran_COMPILER_ID}|${CMAKE_C_COMPILER_ID}" STREQUAL "Intel|Intel")
    set(_toolchain "intel")
  elseif("${CMAKE_Fortran_COMPILER_ID}|${CMAKE_C_COMPILER_ID}" STREQUAL "NAG|GNU")
    set(_toolchain "nag")
  else()
    set(_toolchain "generic")
  endif()
    
  set(${toolchain} "${_toolchain}" PARENT_SCOPE)
  
endfunction()


# Loads toolchain settings.
#
macro(dftbp_load_toolchain_settings)
  
  if(NOT DEFINED TOOLCHAIN_FILE AND NOT "$ENV{DFTBPLUS_TOOLCHAIN_FILE}" STREQUAL "")
    set(TOOLCHAIN_FILE "$ENV{DFTBPLUS_TOOLCHAIN_FILE}")
  endif()
  if(NOT DEFINED TOOLCHAIN AND NOT "$ENV{DFTBPLUS_TOOLCHAIN}" STREQUAL "")
    set(TOOLCHAIN "$ENV{DFTBPLUS_TOOLCHAIN}")
  endif()
  if(NOT DEFINED TOOLCHAIN_FILE OR TOOLCHAIN_FILE STREQUAL "")
    if(NOT DEFINED TOOLCHAIN OR TOOLCHAIN STREQUAL "")
      dftbp_guess_toolchain(TOOLCHAIN)
    endif()
    set(TOOLCHAIN_FILE ${CMAKE_SOURCE_DIR}/sys/${TOOLCHAIN}.cmake)
  endif()
  message(STATUS "Reading build environment specific toolchain file: ${TOOLCHAIN_FILE}")
  include(${TOOLCHAIN_FILE})
endmacro()


# Sets up the global compiler flags
#
macro (dftbp_setup_global_compiler_flags)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILDTYPE_UPPER)
  foreach (lang IN ITEMS Fortran C)
    set(CMAKE_${lang}_FLAGS " ${${lang}_FLAGS}")
    set(CMAKE_${lang}_FLAGS_${BUILDTYPE_UPPER} " ${${lang}_FLAGS_${BUILDTYPE_UPPER}}")
    message(STATUS "Flags for ${lang}-compiler: "
      "${CMAKE_${lang}_FLAGS} ${CMAKE_${lang}_FLAGS_${BUILDTYPE_UPPER}}")
  endforeach()

endmacro()
