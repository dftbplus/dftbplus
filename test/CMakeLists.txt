add_subdirectory(src/dftbp)
add_subdirectory(app)
add_subdirectory(tools/dptools)

#
# General options for all targets
#
set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${CMAKE_CURRENT_SOURCE_DIR}/include -DRELEASE="'${RELEASE}'")

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)
