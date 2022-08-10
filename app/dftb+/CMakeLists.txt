if(NOT BUILD_EXPORTED_TARGETS_ONLY)

  set(projectdir ${PROJECT_SOURCE_DIR})

  set(fypp_flags ${FYPP_BUILD_FLAGS})
  list(APPEND fypp_flags -I${projectdir}/src/dftbp/include -DRELEASE="'${RELEASE}'")

  set(ALL-SOURCES-F90)
  set(ALL-SOURCES-FPP
    dftbplus.F90)

  dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${ALL-SOURCES-FPP}"
    all-sources-f90-preproc)

  add_executable(dftb+ ${ALL-SOURCES-F90} ${all-sources-f90-preproc})

  if (EXTERNAL_MODEL_LIBRARIES)
    target_link_libraries(dftb+ dftbplus ${EXTERNAL_MODEL_LIBRARIES})
  else()
    target_link_libraries(dftb+ dftbplus)
  endif()


  install(TARGETS dftb+ DESTINATION ${CMAKE_INSTALL_BINDIR})

endif()
