set(projectdir ${PROJECT_SOURCE_DIR})

# Collect commonly used modules into object library to prevent .mod file renaming conflicts during
# parallel compilation (see e.g. https://github.com/Unidata/netcdf-fortran/issues/64)

add_library(test_api_mm_utils OBJECT testhelpers.f90 extchargepot.f90 extchargepotgen.f90)
target_link_libraries(test_api_mm_utils PUBLIC dftbplus)

add_executable(test_fileinit test_fileinit.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_fileinit dftbplus)

if(WITH_C_EXECUTABLES)

  add_executable(test_fileinitc test_fileinitc.c $<TARGET_OBJECTS:test_api_mm_utils>)
  target_link_libraries(test_fileinitc dftbplus)
  set_target_properties(test_fileinitc PROPERTIES LINKER_LANGUAGE C)

  add_executable(test_neighbour_list_c test_neighbour_list_c.c $<TARGET_OBJECTS:test_api_mm_utils>)
  target_link_libraries(test_neighbour_list_c dftbplus)
  set_target_properties(test_neighbour_list_c PROPERTIES LINKER_LANGUAGE C)

  add_executable(test_qdepextpotc test_qdepextpotc.c $<TARGET_OBJECTS:test_api_mm_utils>)
  target_link_libraries(test_qdepextpotc dftbplus)
  set_target_properties(test_qdepextpotc PROPERTIES LINKER_LANGUAGE C)

  if(WITH_MPI)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    add_executable(test_mpic test_mpic.c $<TARGET_OBJECTS:test_api_mm_utils>)
    target_link_libraries(test_mpic dftbplus ${MPI_C_LIBRARIES})
    set_target_properties(test_mpic PROPERTIES LINKER_LANGUAGE C)
  endif()

endif()

add_executable(test_treeinit test_treeinit.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_treeinit dftbplus)

add_executable(test_extpot test_extpot.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_extpot dftbplus)

add_executable(test_extpot2 test_extpot2.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_extpot2 dftbplus)

add_executable(test_extcharges test_extcharges.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_extcharges dftbplus)

add_executable(test_extefield test_extefield.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_extefield dftbplus)

add_executable(test_intcharges test_intcharges.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_intcharges dftbplus)

add_executable(test_qdepextpot test_qdepextpot.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_qdepextpot dftbplus)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/src/dftbp/include -DRELEASE="'${releasename}'")

add_executable(test_timeprop test_timeprop.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_timeprop dftbplus)

add_executable(test_ehrenfest test_ehrenfest.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_ehrenfest dftbplus)

add_executable(test_ehrenfest_ext_ions test_ehrenfest_ext_ions.f90
  $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_ehrenfest_ext_ions dftbplus)

add_executable(test_neighbour_list test_neighbour_list.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_neighbour_list dftbplus)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "test_setspeciesanddependents.F90"
  test_setspeciesanddependents.f90)

add_executable(test_setspeciesanddependents test_setspeciesanddependents.f90
  $<TARGET_OBJECTS:test_api_mm_utils>)
target_link_libraries(test_setspeciesanddependents PRIVATE dftbplus)

if(WITH_MPI)
  target_include_directories(test_setspeciesanddependents PRIVATE ${MPI_Fortran_MODULE_DIR})
  target_link_libraries(test_setspeciesanddependents PRIVATE ${MPI_Fortran_LIBRARIES})

  add_executable(test_mpisubgrids test_mpisubgrids.f90 $<TARGET_OBJECTS:test_api_mm_utils>)
  target_link_libraries(test_mpisubgrids PRIVATE dftbplus)
  target_include_directories(test_mpisubgrids PRIVATE ${MPI_Fortran_MODULE_DIR})
  target_link_libraries(test_mpisubgrids PRIVATE ${MPI_Fortran_LIBRARIES})

  set(mpi_sources
    test_setspeciesanddependents.f90
    test_mpisubgrids.f90)

  # Lift strict actual argument type checking to enable mpi-calls to compile
  foreach(mpisource IN LISTS mpi_sources)
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
      set_source_files_properties(SOURCE ${mpisource} PROPERTY COMPILE_FLAGS -mismatch)
    endif()
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
      if("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_GREATER_EQUAL "11")
        set_source_files_properties(SOURCE ${mpisource} PROPERTY COMPILE_FLAGS
          -fallow-argument-mismatch)
      endif()
    endif()
  endforeach()

endif()
