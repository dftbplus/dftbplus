set(sources-fpp
  initmodes.F90
  modeprojection.F90
  modes.F90)

list(APPEND FYPP_FLAGS -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

dftbp_preprocess("${FYPP}" "${FYPP_FLAGS}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(modes ${sources-f90-preproc})

if(WITH_MPI)
  target_include_directories(modes PRIVATE ${MPI_Fortran_MODULE_DIR})
  target_link_libraries(modes PRIVATE ${MPI_Fortran_LIBRARIES})
endif()

target_link_libraries(modes PRIVATE dftbplus)


install(TARGETS modes DESTINATION ${INSTALL_BIN_DIR})
