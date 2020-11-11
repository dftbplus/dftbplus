file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../" projectdir)

set(sources-fpp
  libnegf_int.F90
  initphonons.F90
  phonons.F90)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/prog/dftb+/include)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(phonons ${sources-f90-preproc})

if(WITH_MPI)
  target_link_libraries(phonons PRIVATE MPI::MPI_Fortran)
endif()

target_link_libraries(phonons PUBLIC dftbplus)
target_link_libraries(phonons PUBLIC negf)

install(TARGETS phonons DESTINATION ${CMAKE_INSTALL_BINDIR})
