file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../../" projectdir)

# Collect commonly used modules into object library to prevent .mod file renaming conflicts during
# parallel compilation (see e.g. https://github.com/Unidata/netcdf-fortran/issues/64)
add_library(test_api_mm_utils OBJECT testhelpers.f90 extchargepot.f90 extchargepotgen.f90)
add_dependencies(test_api_mm_utils dftbplus)
target_include_directories(test_api_mm_utils PUBLIC
  $<TARGET_PROPERTY:dftbplus,INTERFACE_INCLUDE_DIRECTORIES>)

add_executable(test_fileinit test_fileinit.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_fileinit dftbplus)

if(WITH_C_EXECUTABLES)
  add_executable(test_fileinitc test_fileinitc.c $<TARGET_OBJECTS:test_api_mm_utils>)
  target_link_libraries(test_fileinitc dftbplus)
  set_target_properties(test_fileinitc PROPERTIES LINKER_LANGUAGE C)
endif()

add_executable(test_treeinit test_treeinit.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_treeinit dftbplus)

add_executable(test_extpot test_extpot.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_extpot dftbplus)

add_executable(test_extpot2 test_extpot2.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_extpot2 dftbplus)

add_executable(test_extcharges test_extcharges.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_extcharges dftbplus)

add_executable(test_qdepextpot test_qdepextpot.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_qdepextpot dftbplus)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/prog/dftb+/include -DRELEASE="'${releasename}'")

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "test_setspeciesanddependents.F90"
  test_setspeciesanddependents.f90)

# NAG compiler won't compile files with mpi calls without the '-mismatch' option
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
  set_source_files_properties(SOURCE test_setspeciesanddependents.f90
    PROPERTY COMPILE_FLAGS -mismatch)
endif()

add_executable(test_setspeciesanddependents test_setspeciesanddependents.f90
  $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_setspeciesanddependents PRIVATE dftbplus)
if(WITH_MPI)
  target_include_directories(test_setspeciesanddependents PRIVATE ${MPI_Fortran_MODULE_DIR})
  target_link_libraries(test_setspeciesanddependents PRIVATE ${MPI_Fortran_LIBRARIES})
endif()


if(WITH_C_EXECUTABLES)
  add_executable(test_qdepextpotc test_qdepextpotc.c $<TARGET_OBJECTS:test_api_mm_utils>)
  target_link_libraries(test_qdepextpotc dftbplus)
  set_target_properties(test_qdepextpotc PROPERTIES LINKER_LANGUAGE C)
endif()
