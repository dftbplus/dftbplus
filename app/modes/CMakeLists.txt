set(projectdir ${PROJECT_SOURCE_DIR})

set(sources-fpp
  hsdhelpers.F90
  initmodes.F90
  inputdata.F90
  main.F90
  modeprojection.F90
  modes.F90
  parser.F90)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(modes ${sources-f90-preproc})

if(WITH_MPI)
  target_link_libraries(modes PRIVATE MPI::MPI_Fortran)
endif()

target_link_libraries(modes PRIVATE dftbplus)

install(TARGETS modes DESTINATION "${CMAKE_INSTALL_BINDIR}")
