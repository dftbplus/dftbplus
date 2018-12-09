set(BUILDDIR ${CMAKE_CURRENT_BINARY_DIR})
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR})

_get_fypp_defines(FYPP-OPTS)
list(APPEND FYPP-OPTS -I${CMAKE_SOURCE_DIR}/prog/dftb+/include -DRELEASE="'${RELEASE}'")

add_custom_command(
  OUTPUT _dftbplus_tests
  COMMAND ${CMAKE_SOURCE_DIR}/utils/test/testlist_to_fypp < ${SRCDIR}/tests > ${BUILDDIR}/_dftbplus_tests.fypp
  COMMAND ${FYPP} ${FYPP-OPTS} -DMPI_PROCS=${TEST_MPI_PROCS} -DOMP_THREADS=${TEST_OMP_THREADS} ${BUILDDIR}/tests.fypp > ${BUILDDIR}/_dftbplus_tests
  BYPRODUCTS _dftbplus_tests.fypp
  DEPENDS ${SRCDIR}/tests)

add_custom_target(dftbplus_tests ALL DEPENDS _dftbplus_tests)

add_test(
  NAME test_dftb+
  COMMAND cmake -E echo "Hello")

set_tests_properties(test_dftb+ PROPERTIES DEPENDS _dftbplus_tests)
