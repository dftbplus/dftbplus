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
function(dftbp_register_preprocessing preproc preprocopts oldext newext oldfiles newfiles)

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


# Registers a target as install target with default install directory locations
#
# Args:
#     target [in]: Target to register for installation
#
function(dftbp_register_install_target target)

  install(TARGETS ${target}
    RUNTIME DESTINATION ${INSTALL_BIN_DIR}
    LIBRARY DESTINATION ${INSTALL_LIB_DIR}
    ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
    INCLUDES DESTINATION ${INSTALL_INC_DIR})

endfunction()


# Registers a directory with modfiles for installation
function(dftbp_register_install_mod_dirs moddirs)

  foreach(moddir IN LISTS moddirs)
    install(
      DIRECTORY "${moddir}/"
      DESTINATION "${INSTALL_MOD_DIR}")
  endforeach()

endfunction()


# Build -D command line arguments for Fypp preprocessor based on current settings
#
# Args:
#     fyppdefines [out]: Command line option with -D defines
#
function (dftbp_get_fypp_defines fyppdefines)

  set(_fyppdefines)

  if(WITH_ARPACK)
    list(APPEND _fyppdefines -DWITH_ARPACK)
  endif()

  if(WITH_DFTD3)
    list(APPEND _fyppdefines -DWITH_DFTD3)
  endif()

  if(WITH_MPI)
    list(APPEND _fyppdefines -DWITH_MPI -DWITH_SCALAPACK)
  endif()

  if(WITH_SOCKETS)
    list(APPEND _fyppdefines -DWITH_SOCKETS)
  endif()

  set(${fyppdefines} ${_fyppdefines} PARENT_SCOPE)

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


# Logs settings to file
#
# Args:
#     file [in]: Name of the file to write the log to.
#     * [in]: Any variable name, for which name and value should be logged
#
function (dftbp_log_settings file)

  foreach(var IN LISTS ARGN)
    set(msg "SET(${var} \"${${var}}\")")
    file(APPEND ${file} "${msg}\n")
  endforeach()

endfunction()


# Resets a file to have empty content
#
# Args:
#     file [in]: Name of the file to reset.
#
function (dftbp_reset_file file)
  file(WRITE ${file} "")
endfunction()
