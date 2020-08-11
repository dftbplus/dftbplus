add_executable(buildwire buildwire.f90)
install(TARGETS buildwire DESTINATION ${INSTALL_BIN_DIR})

add_executable(flux flux.f90)
install(TARGETS flux DESTINATION ${INSTALL_BIN_DIR})

add_executable(makecube makecube.f90)
install(TARGETS makecube DESTINATION ${INSTALL_BIN_DIR})

# Sources to be pre-processed for setupgeom
set(setupgeom-fpp
  helpsetupgeom.F90
  inputdata.F90
  parser_setup.F90
  setupgeom.F90)

list(APPEND FYPP_FLAGS -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

dftbp_preprocess("${FYPP}" "${FYPP_FLAGS}" "F90" "f90" "${setupgeom-fpp}" sources-f90-preproc)

add_executable(setupgeom ${sources-f90-preproc})
target_link_libraries(setupgeom dftbplus)

install(TARGETS setupgeom DESTINATION ${INSTALL_BIN_DIR})
