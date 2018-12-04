set(SOURCES-FPP
  gridcache.F90
  initwaveplot.F90
  molorb.F90
  slater.F90
  waveplot.F90)

set(FYPP-OPTS -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

_register_custom_preprocessing(
  "${FYPP}" "${FYPP-OPTS}" "F90" "f90" "${SOURCES-FPP}" SOURCES-F90-PREPROC)

add_executable("waveplot" ${SOURCES-F90-PREPROC})

target_link_libraries(waveplot dftbplus ${LAPACK_LIBRARIES})

_register_install_target(waveplot)
