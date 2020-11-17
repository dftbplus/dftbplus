# Finds customized libraries.
#
# The customized list of libraries can contain library names, file names or linker options.
# Libraries are searched for by the find_library() function, file names are checked on existence and
# linker options (names starting with "-") are added unchanged.
#
# libs [in]: List of libraries to find.
# libdirs [in]: Directories to look up for the libraries.
# find_quietly [in]: Whether the libraries should be found quietly.
# libs_found [out]: Contains the list of the detected libraries with their full path. If any of
#     the libraries/files could not be found, the corresponding entry is replaced by
#     "_libpath-NOT_FOUND".
#
function(find_custom_libraries libs libdirs find_quietly libs_found)

  set(_libs)
  foreach(_lib IN LISTS libs)
    if(_lib MATCHES "^[ ]*-.*"  OR EXISTS "${_lib}")
      list(APPEND _libs ${_lib})
    else()
      find_library(_libpath ${_lib} HINTS ${libdirs})
      list(APPEND _libs ${_libpath})
      if(NOT _libpath AND NOT find_quietly)
        message(WARNING "Could not find library '${_lib}' (path hints '${libdirs}'")
      endif()
      unset(_libpath CACHE)
    endif()
  endforeach()
  set(${libs_found} "${_libs}" PARENT_SCOPE)

endfunction()
