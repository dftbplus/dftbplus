set(projectdir ${PROJECT_SOURCE_DIR})

set(sources-fpp printunits.F90)

set(fypp-opts -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(printunits ${sources-f90-preproc})

target_link_libraries(printunits PRIVATE dftbplus)

install(TARGETS printunits DESTINATION "${CMAKE_INSTALL_BINDIR}")
