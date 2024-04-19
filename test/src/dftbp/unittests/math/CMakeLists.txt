set(unit-test-prefix "${unit-test-prefix}/math")

if (NOT WITH_MPI)
#  set(tests
#    matrixoperations_mpi)
#  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
#else()
  set(tests
    matrixoperations)
endif()

foreach(test IN LISTS tests)
  fytest_fypp_preprocess("${FYPP_FLAGS}" ${test}.F90 sources-f90)
  add_executable(${test} ${sources-f90})
  target_link_libraries(${test} dftbplus FyTest::FyTest)
  #if (WITH_MPI)
  #  target_include_directories(${test} PRIVATE ${MPI_Fortran_MODULE_DIR})
  #  target_link_libraries(${test} MPI::MPI_Fortran MpiFx::MpiFx
  #    ScalapackFx::ScalapackFx ${MPI_Fortran_LIBRARIES})
  #endif()
  add_test(NAME ${unit-test-prefix}/${test} COMMAND ${test})
endforeach()
