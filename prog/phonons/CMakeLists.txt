set(sources-fpp
  libnegf_int.F90
  initprogram.F90
  phonons.F90)

set(fypp-opts -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

dftbp_preprocess("${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(phonons ${sources-f90-preproc})
target_link_libraries(phonons dftbplus)

install(TARGETS phonons DESTINATION ${INSTALL_BIN_DIR})
