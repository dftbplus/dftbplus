set(builddir ${CMAKE_CURRENT_BINARY_DIR})
set(srcdir ${CMAKE_CURRENT_SOURCE_DIR})
file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../" projectdir)

# execute_process is verbatim, therefore ";" must be replaced with " "
string(REGEX REPLACE ";" " " releasename "${RELEASE}")


execute_process(
  COMMAND ${projectdir}/utils/test/testlist_to_fypp
  INPUT_FILE ${srcdir}/tests
  OUTPUT_FILE ${builddir}/_dftbplus_tests.fypp)

set(fypp_flags ${FYPP_CONFIG_FLAGS})
list(APPEND fypp_flags -I${projectdir}/prog/dftb+/include -DRELEASE="'${releasename}'")

execute_process(
  COMMAND ${FYPP} ${fypp_flags} -DMPI_PROCS=${TEST_MPI_PROCS} -DOMP_THREADS=${TEST_OMP_THREADS}
  INPUT_FILE ${builddir}/_dftbplus_tests.fypp
  OUTPUT_FILE ${builddir}/_dftbplus_tests)

file(STRINGS ${builddir}/_dftbplus_tests dftbplus_tests_raw)
foreach(line IN LISTS dftbplus_tests_raw)
  string(STRIP "${line}" testname)
  if(NOT "${testname}" STREQUAL "")
    list(APPEND tests ${testname})
  endif()
endforeach()

foreach(test IN LISTS tests)
  add_test(
    NAME dftb+_${test}
    COMMAND ${srcdir}/bin/autotest2
        -r ${srcdir} -w ${builddir} -d ${srcdir}/bin/tagdiff
        -P "${TEST_RUNNER}" -p "$<TARGET_FILE:dftb+>"
        -s P,R,C,S ${test})
endforeach()
