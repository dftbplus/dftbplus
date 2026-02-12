set(builddir ${CMAKE_CURRENT_BINARY_DIR})
set(srcdir ${CMAKE_CURRENT_SOURCE_DIR})
set(projectdir ${PROJECT_SOURCE_DIR})

execute_process(
  COMMAND ${projectdir}/utils/test/testlist_to_fypp
  INPUT_FILE ${srcdir}/tests
  OUTPUT_FILE ${builddir}/_dftbplus_tests.fypp)

set(fypp_flags ${FYPP_CONFIG_FLAGS})
list(APPEND fypp_flags -I${projectdir}/src/dftbp/include)

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
  if(EXISTS ${srcdir}/${test}/testecute.toml)
    set(testecute_toml "${srcdir}/${test}/testecute.toml")
  else()
    set(testecute_toml "${srcdir}/testecute.toml")
  endif()
  set(tagdiffconfopts "-c ${srcdir}/tagdiff.toml")
  if(EXISTS ${srcdir}/${test}/tagdiff.toml)
    set(tagdiffconfopts "-c ${srcdir}/${test}/tagdiff.toml ${tagdiffconfopts}")
  endif()
  file(MAKE_DIRECTORY ${builddir}/${test})
  add_test(
    NAME app/dftb+/${test}
    COMMAND
        ${PYTHON_INTERPRETER} -m testecute
        -P ${srcdir}/testecute_plugins
        -D testdir=${srcdir}/${test}
        -D workdir=${builddir}/${test}
        -D executable=${TEST_RUNNER}\ $<TARGET_FILE:dftb+>
        -D paramdir=${PROJECT_SOURCE_DIR}/external
        -D tagdiff=${srcdir}/bin/tagdiff
        -D tagdiffconfopts=${tagdiffconfopts}
        ${testecute_toml}
  )
  set_tests_properties(
    app/dftb+/${test}
    PROPERTIES
    WORKING_DIRECTORY
      ${builddir}/${test}
    ENVIRONMENT
      "PYTHONPATH=${PROJECT_SOURCE_DIR}/external/testecute/origin/src")
endforeach()
