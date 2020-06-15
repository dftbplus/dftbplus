set(sources-fpp
  helpsetupgeom.F90
  inputdata.F90
  parser_setup.F90
  setupgeom.F90)

list(APPEND FYPP_FLAGS -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

dftbp_preprocess("${FYPP}" "${FYPP_FLAGS}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(setupgeom ${sources-f90-preproc})
target_link_libraries(setupgeom dftbplus)

install(TARGETS setupgeom DESTINATION ${INSTALL_BIN_DIR})
