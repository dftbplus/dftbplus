set(SOURCES-FPP
  initmodes.F90
  modeprojection.F90
  modes.F90)

set(FYPP-OPTS -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

_register_custom_preprocessing(
  "${FYPP}" "${FYPP-OPTS}" "F90" "f90" "${SOURCES-FPP}" SOURCES-F90-PREPROC)

add_executable("modes" ${SOURCES-F90-PREPROC})

target_link_libraries(modes dftbplus ${LAPACK_LIBRARIES})

_register_install_target(modes)
