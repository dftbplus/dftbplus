set(unit-test-prefix "${unit-test-prefix}/io")

set(tests
  indexselection
  tokenreader)

foreach(test IN LISTS tests)
  fytest_fypp_preprocess("${FYPP_FLAGS}" ${test}.F90 sources-f90)
  add_executable(${test} ${sources-f90})
  target_link_libraries(${test} dftbplus FyTest::FyTest)
  add_test(NAME ${unit-test-prefix}/${test} COMMAND ${test})
endforeach()
