set(sources-fpp skderivs.F90)

set(fypp-opts -I${CMAKE_SOURCE_DIR}/prog/dftb+/include)

dftbp_preprocess("${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(skderivs ${sources-f90-preproc})

target_link_libraries(skderivs PRIVATE dftbplus)

install(TARGETS skderivs DESTINATION ${INSTALL_BIN_DIR})
