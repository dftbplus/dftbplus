set(unit-test-prefix "${unit-test-prefix}/include")

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${PROJECT_SOURCE_DIR}/src/dftbp/include)

set(test allocatablelist)
set(sources allocatablelist.F90 alloclisthelper.F90)
fytest_fypp_preprocess("${fypp_flags}" "${sources}" sources-f90)
add_executable(${test} ${sources-f90})
target_link_libraries(${test} dftbplus FyTest::FyTest)
add_test(NAME ${unit-test-prefix}/${test} COMMAND ${test})

set(test pointerlist)
set(sources pointerlist.F90 ptrlisthelper.F90)
fytest_fypp_preprocess("${fypp_flags}" "${sources}" sources-f90)
add_executable(${test} ${sources-f90})
target_link_libraries(${test} dftbplus FyTest::FyTest)
add_test(NAME ${unit-test-prefix}/${test} COMMAND ${test})
