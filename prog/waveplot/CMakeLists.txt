set(sources-fpp
  gridcache.F90
  initwaveplot.F90
  molorb.F90
  slater.F90
  waveplot.F90)

set(fypp-opts -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

_register_custom_preprocessing(
  "${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable("waveplot" ${sources-f90-preproc})

target_link_libraries(waveplot dftbplus ${LAPACK_LIBRARIES})

_register_install_target(waveplot)
