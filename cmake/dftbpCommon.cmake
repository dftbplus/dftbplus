# Replaces the extension of a given file
function(_replace_extension oldext newext fname newfname)

  string(REGEX REPLACE "\\.${oldext}$" ".${newext}" _newfname ${fname})
  set(${newfname} ${_newfname} PARENT_SCOPE)

endfunction(_replace_extension)


# Registers files for preprocessing
function(_register_custom_preprocessing prep prepopts oldext newext oldfiles newfiles)

  set(_newfiles)
  foreach(oldfile ${oldfiles})
    _replace_extension(${oldext} ${newext} ${oldfile} newfile)
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${newfile}
      COMMAND ${prep} ${prepopts} ${CMAKE_CURRENT_SOURCE_DIR}/${oldfile} ${CMAKE_CURRENT_BINARY_DIR}/${newfile}
      MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${oldfile})
    list(APPEND _newfiles ${CMAKE_CURRENT_BINARY_DIR}/${newfile})
  endforeach()
  set(${newfiles} ${_newfiles} PARENT_SCOPE)

endfunction(_register_custom_preprocessing)


# Registers a target as to be installed.
function(_register_install_target target)

  install(TARGETS ${target}
    RUNTIME DESTINATION ${INSTALL_BIN_DIR}
    LIBRARY DESTINATION ${INSTALL_LIB_DIR}
    ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
    INCLUDES DESTINATION ${INSTALL_INC_DIR})

endfunction(_register_install_target)


# Registers a directory with modfiles for installation
function(_register_install_modfile_dirs moddirs)

  foreach(moddir ${moddirs})
    install(
      DIRECTORY "${moddir}/"
      DESTINATION "${INSTALL_MOD_DIR}")
  endforeach()

endfunction(_register_install_modfile_dirs)


# Build Fypp -D command line arguments based on current settings
function (_get_fypp_defines fyppdefines)

  set(_fyppdefines)

  if(WITH_ARPACK)
    list(APPEND _fyppdefines -DWITH_ARPACK)
  endif()

  if(WITH_DFTD3)
    list(APPEND _fyppdefines -DWITH_DFTD3)
  endif()

  if(WITH_MPI)
    list(APPEND _fyppdefines -DWITH_MPI -DWITH_SCALAPACK)
  endif(WITH_MPI)

  if(WITH_SOCKETS)
    list(APPEND _fyppdefines -DWITH_SOCKETS)
  endif(WITH_SOCKETS)

  set(${fyppdefines} ${_fyppdefines} PARENT_SCOPE)

endfunction(_get_fypp_defines)


# Gets the release information
function(_get_release_name release)

  if (NOT EXISTS ${CMAKE_BINARY_DIR}/RELEASE)
    if(EXISTS ${CMAKE_SOURCE_DIR}/RELEASE)
      file(COPY ${CMAKE_SOURCE_DIR}/RELEASE DESTINATION ${CMAKE_BINARY_DIR})
    else()
      execute_process(
	COMMAND ${CMAKE_SOURCE_DIR}/utils/build/update_release ${CMAKE_BINARY_DIR}/RELEASE
	RESULT_VARIABLE exitcode)
      if(NOT exitcode EQUAL "0")
	file(WRITE ${CMAKE_BINARY_DIR}/RELEASE "(UNKNOWN RELEASE)")
      endif()
    endif()
  endif()
  file(READ ${CMAKE_BINARY_DIR}/RELEASE _release)
  separate_arguments(_release)
  set(${release} "${_release}" PARENT_SCOPE)

endfunction(_get_release_name)


# Logs settings
function (_log_settings file variables)

  foreach (var ${variables})
    set(msg "SET(${var} ${${var}})")
    string(REPLACE ";" " " msgstr "${msg}")
    file(APPEND ${file} "${msgstr}\n")
  endforeach()

endfunction(_log_settings)


function (_reset_file file)
  file(WRITE ${file} "")
endfunction (_reset_file)
