set(targets
  integvalue
  polyvalue
  splvalue)

set(fypp-opts -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

foreach(target ${targets})

  set(sources-fpp "${target}.F90")

  _register_custom_preprocessing(
    "${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

  add_executable(${target} ${sources-f90-preproc})

  target_link_libraries("${target}" dftbplus ${LAPACK_LIBRARIES})

  _register_install_target(${target})

endforeach()
