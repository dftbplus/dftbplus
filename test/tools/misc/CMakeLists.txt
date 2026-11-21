if (WITH_PYTHON)
  add_test(
    NAME test_mo_irreps
    COMMAND ${PYTHON_INTERPRETER} ${CMAKE_CURRENT_SOURCE_DIR}/test_mo_irreps.py)
endif()
