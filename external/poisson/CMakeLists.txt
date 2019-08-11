set(sources-fpp
  bulkpot.F90
  fancybc.F90
  gallocation.F90
  gclock.F90
  gconstants.F90
  gewald.F90
  gprecision.F90
  mpi_poisson.F90
  parameters.F90
  parcheck.F90
  poisson.F90
  std_io.F90
  structure.F90)

set(sources-f77
  mudpack/dmud3.f
  mudpack/dmuh3.f
  mudpack/dmud3ln.f
  mudpack/dmud3pn.f
  mudpack/dmud3sp.f
  mudpack/dmudcom.f)

dftbp_register_preprocessing(
  "${FYPP}" "${FYPP_FLAGS}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_library(poisson ${sources-f90-preproc} ${sources-f77})

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(poisson PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})
target_include_directories(poisson PUBLIC ${includedir})

if(WITH_MPI)
  target_link_libraries(poisson PUBLIC mpifx)
endif()
