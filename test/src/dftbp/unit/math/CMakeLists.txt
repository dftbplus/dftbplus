set(test-collection "math")

if (NOT WITH_MPI)

  set(sources-f90
    testapp.f90)

  set(sources-fypp
    interpolation.F90
    matrixops.F90)

  set(sources-fypp-f90)
  foreach(source-fypp IN LISTS sources-fypp)
    dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${source-fypp}" source-fypp-f90)
    list(APPEND sources-fypp-f90 ${source-fypp-f90})
  endforeach()

  set(testapp "testapp_${test-collection}")
  add_executable(${testapp})
  target_sources(${testapp} PRIVATE ${sources-f90} ${sources-fypp-f90})
  target_link_libraries(${testapp} dftbplus Fortuno::fortuno_serial)
  add_test(NAME ${unit-test-prefix}/${test-collection} COMMAND ${testapp})

endif ()
