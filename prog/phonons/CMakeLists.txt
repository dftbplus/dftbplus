set(sources-fpp
  masses.F90	
  libnegf_int.F90
  initphonons.F90
  phonons.F90)

set(fypp-opts -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

dftbp_preprocess("${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(phonons ${sources-f90-preproc})
target_link_libraries(phonons PUBLIC dftbplus)
target_link_libraries(phonons PUBLIC negf)

install(TARGETS phonons DESTINATION ${CMAKE_INSTALL_BINDIR})
