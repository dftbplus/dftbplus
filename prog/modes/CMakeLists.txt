set(sources-fpp
  initmodes.F90
  modeprojection.F90
  modes.F90)

set(fypp-opts -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

dftbp_preprocess("${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(modes ${sources-f90-preproc})

target_link_libraries(modes dftbplus)

install(TARGETS modes DESTINATION ${INSTALL_BIN_DIR})
