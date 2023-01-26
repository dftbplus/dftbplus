if (WITH_PYTHON)
  add_test(
    NAME test_ptools
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/runtests.sh ${PYTHON_INTERPRETER})
  set_tests_properties(
    test_ptools
    PROPERTIES
    ENVIRONMENT "DFTBPLUS_PARAM_DIR=${PROJECT_SOURCE_DIR}/external")
endif()
