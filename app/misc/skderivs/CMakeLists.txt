set(projectdir ${PROJECT_SOURCE_DIR})

set(sources-fpp skderivs.F90)

set(fypp-opts -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(skderivs ${sources-f90-preproc})

target_link_libraries(skderivs PRIVATE dftbplus)

install(TARGETS skderivs DESTINATION "${CMAKE_INSTALL_BINDIR}")
