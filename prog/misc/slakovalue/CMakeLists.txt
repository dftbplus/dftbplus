set(targets
  integvalue
  polyvalue
  splvalue)

set(fypp-opts -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

foreach(target IN LISTS targets)

  set(sources-fpp "${target}.F90")

  dftbp_preprocess("${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

  add_executable(${target} ${sources-f90-preproc})

  target_link_libraries(${target} dftbplus)

  install(TARGETS ${target} DESTINATION ${INSTALL_BIN_DIR})

endforeach()
