set(unit-test-prefix "${unit-test-prefix}/include")

set(tests
  pointerlist
  allocatablelist)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${PROJECT_SOURCE_DIR}/src/dftbp/include -DUNIT_TEST=True)

foreach(test IN LISTS tests)
  fytest_fypp_preprocess("${fypp_flags}" ${test}.F90 sources-f90)
  add_executable(${test} ${sources-f90})
  target_link_libraries(${test} dftbplus FyTest::FyTest)
  add_test(NAME ${unit-test-prefix}/${test} COMMAND ${test})
endforeach()
