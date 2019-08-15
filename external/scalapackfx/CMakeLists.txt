# Make sure, we do not overwrite something which is relevant for the host project
if(DEFINED LIBRARY_ONLY)
  message(FATAL_ERROR "Possible conflict of cache variable LIBRARY_ONLY")
endif()

set(LIBRARY_ONLY TRUE CACHE BOOL "Whether only the library should be compiled")

if(BUILD_API)
  add_subdirectory(origin)
else()
  add_subdirectory(origin EXCLUDE_FROM_ALL)
endif()

unset(LIBRARY_ONLY CACHE)
