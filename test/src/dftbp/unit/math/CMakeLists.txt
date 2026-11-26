set(test-collection "math")

if (NOT WITH_MPI)

  set(sources-fypp
    binarysearch.F90
    lapacksvd.F90
    matrixops.F90
    sorting.F90)

  set(sources-f90
    testapp.f90)

else()

  set(sources-fypp
    mpi_matrixops.F90)

  set(sources-f90
    mpi_testapp.f90)

endif ()

set(sources-fypp-f90)
foreach(source-fypp IN LISTS sources-fypp)
  dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${source-fypp}" source-fypp-f90)
  list(APPEND sources-fypp-f90 ${source-fypp-f90})
endforeach()

set(testapp "testapp_${test-collection}")
add_executable(${testapp})
target_sources(${testapp} PRIVATE ${sources-f90} ${sources-fypp-f90})

if (WITH_MPI)

  target_link_libraries(${testapp} dftbplus Fortuno::fortuno_mpi MPI::MPI_Fortran)

  add_test(
    NAME ${unit-test-prefix}/${test-collection}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_MPI_PROCS}
            ${MPIEXEC_PREFLAGS} ${testapp} ${MPIEXEC_POSTFLAGS}
  )


else()

  target_link_libraries(${testapp} dftbplus Fortuno::fortuno_serial)
  add_test(NAME ${unit-test-prefix}/${test-collection} COMMAND ${testapp})

endif ()
