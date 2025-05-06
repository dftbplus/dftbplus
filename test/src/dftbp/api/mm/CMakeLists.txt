set(projectdir ${PROJECT_SOURCE_DIR})
set(builddir ${CMAKE_CURRENT_BINARY_DIR})
set(srcdir ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(testers)

set(fypp_flags ${FYPP_CONFIG_FLAGS} -I${projectdir}/src/dftbp/include)

execute_process(
  COMMAND ${projectdir}/utils/test/testlist_to_fypp
  INPUT_FILE ${srcdir}/tests
  OUTPUT_FILE ${builddir}/_api_tests.fypp)

execute_process(
  COMMAND ${FYPP} ${fypp_flags} -DMPI_PROCS=${TEST_MPI_PROCS} -DOMP_THREADS=${TEST_OMP_THREADS}
  INPUT_FILE ${builddir}/_api_tests.fypp
  OUTPUT_FILE ${builddir}/_api_tests)

file(STRINGS ${builddir}/_api_tests api_tests_raw)
foreach(testname_raw IN LISTS api_tests_raw)
  string(STRIP "${testname_raw}" testname)
  if(NOT "${testname}" STREQUAL "")
    list(APPEND tests ${testname})
  endif()
endforeach()

foreach(test IN LISTS tests)
  add_test(
    NAME api_${test}
    COMMAND ${projectdir}/test/app/dftb+/bin/autotest2
        -r ${srcdir}/testcases -w ${builddir}/testcases
        -d ${projectdir}/test/app/dftb+/bin/tagdiff
        -P "${TEST_RUNNER}" -s P,R,C,S ${test})
  set_tests_properties(
    api_${test}
    PROPERTIES
    ENVIRONMENT "DFTBPLUS_PARAM_DIR=${PROJECT_SOURCE_DIR}/external")
endforeach()
