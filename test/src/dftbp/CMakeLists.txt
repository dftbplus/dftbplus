if(WITH_UNIT_TESTS AND NOT WITH_MPI)
  add_subdirectory(unit)
endif()

if(WITH_API)
  add_subdirectory(api)
endif()
