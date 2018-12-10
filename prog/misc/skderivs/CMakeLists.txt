set(sources-fpp skderivs.F90)

set(fypp-opts -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

_register_custom_preprocessing(
  "${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable("skderivs" ${sources-f90-preproc})

target_link_libraries(skderivs dftbplus ${LAPACK_LIBRARIES})

_register_install_target(skderivs)
