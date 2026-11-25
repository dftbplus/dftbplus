if(WITH_UNIT_TESTS)
  add_subdirectory(unit)
endif()

if(WITH_API)
  add_subdirectory(api)
endif()
