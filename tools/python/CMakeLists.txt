add_subdirectory(ptools)
if(WITH_API AND ENABLE_DYNAMIC_LOADING)
  add_subdirectory(api)
endif()
