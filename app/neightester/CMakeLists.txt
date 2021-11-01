set(projectdir ${PROJECT_SOURCE_DIR})

set(sources-fpp
    test_neighlist.F90)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(test_neighlist ${sources-f90-preproc})
target_link_libraries(test_neighlist PRIVATE dftbplus)
