if (WITH_PYTHON)
  add_test(
    NAME test_ptools
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtests.sh ${PYTHON_INTERPRETER})
endif()
