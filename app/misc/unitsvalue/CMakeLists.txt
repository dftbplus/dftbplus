set(projectdir ${PROJECT_SOURCE_DIR})

set(sources-fpp print_units.F90)

set(fypp-opts -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp-opts}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(print_units ${sources-f90-preproc})

target_link_libraries(print_units PRIVATE dftbplus)

install(TARGETS print_units DESTINATION "${CMAKE_INSTALL_BINDIR}")
