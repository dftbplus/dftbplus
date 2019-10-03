set(builddir ${CMAKE_CURRENT_BINARY_DIR})
set(srcdir ${CMAKE_CURRENT_SOURCE_DIR})

# execute_process is verbatim, therefore ";" must be replaced with " "
string(REGEX REPLACE ";" " " releasename "${RELEASE}")

list(APPEND FYPP_FLAGS -I${CMAKE_SOURCE_DIR}/prog/dftb+/include -DRELEASE="'${releasename}'")

add_subdirectory(testers)

execute_process(
  COMMAND ${CMAKE_SOURCE_DIR}/utils/test/testlist_to_fypp
  INPUT_FILE ${srcdir}/tests
  OUTPUT_FILE ${builddir}/_api_tests.fypp)

execute_process(
  COMMAND ${FYPP} ${FYPP_FLAGS} -DMPI_PROCS=${TEST_MPI_PROCS} -DOMP_THREADS=${TEST_OMP_THREADS}
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
    COMMAND ${CMAKE_SOURCE_DIR}/test/prog/dftb+/bin/autotest2
        -r ${srcdir}/testcases -w ${builddir}/testcases
	-d ${CMAKE_SOURCE_DIR}/test/prog/dftb+/bin/tagdiff
        -P "${TEST_RUNNER}" -s P,R,C,S ${test})
endforeach()
