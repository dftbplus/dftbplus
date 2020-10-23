file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../" projectdir)

set(sources-fpp
  gridcache.F90
  initwaveplot.F90
  molorb.F90
  slater.F90
  waveplot.F90)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/prog/dftb+/include)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(waveplot ${sources-f90-preproc})

if(WITH_MPI)
  target_link_libraries(waveplot PRIVATE MPI::MPI_Fortran)
endif()

target_link_libraries(waveplot PRIVATE dftbplus)

install(TARGETS waveplot DESTINATION "${CMAKE_INSTALL_BINDIR}")
