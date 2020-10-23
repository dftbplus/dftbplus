file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../" projectdir)

add_executable(buildwire buildwire.f90)
install(TARGETS buildwire DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(flux flux.f90)
install(TARGETS flux DESTINATION "${CMAKE_INSTALL_BINDIR}")

add_executable(makecube makecube.f90)
install(TARGETS makecube DESTINATION "${CMAKE_INSTALL_BINDIR}")

# Sources to be pre-processed for setupgeom
set(setupgeom-fpp
  helpsetupgeom.F90
  inputdata.F90
  parser_setup.F90
  setupgeom.F90)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/prog/dftb+/include)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${setupgeom-fpp}" sources-f90-preproc)

add_executable(setupgeom ${sources-f90-preproc})

if(WITH_MPI)
  target_link_libraries(setupgeom PRIVATE MPI::MPI_Fortran)
endif()

target_link_libraries(setupgeom PRIVATE dftbplus)

install(TARGETS setupgeom DESTINATION "${CMAKE_INSTALL_BINDIR}")
