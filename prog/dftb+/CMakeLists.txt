#
# General options for all targets
#

_get_fypp_defines(FYPP-OPTS)
list(APPEND FYPP-OPTS -I${CMAKE_CURRENT_SOURCE_DIR}/include -DRELEASE="'${RELEASE}'")


#
# Compile and install library
#

set(ALL-SOURCES-F90)
set(ALL-SOURCES-FPP)

add_subdirectory(lib_common)
add_subdirectory(lib_derivs)
add_subdirectory(lib_dftb)
add_subdirectory(lib_extlibs)
add_subdirectory(lib_geoopt)
add_subdirectory(lib_io)
add_subdirectory(lib_math)
add_subdirectory(lib_md)
add_subdirectory(lib_mixer)
add_subdirectory(lib_timedep)
add_subdirectory(lib_type)

_register_custom_preprocessing(
  "${FYPP}" "${FYPP-OPTS}" "F90" "f90" "${ALL-SOURCES-FPP}" ALL-SOURCES-F90-PREPROC)

add_library(dftbplus
  ${ALL-SOURCES-F90}
  ${ALL-SOURCES-F90-PREPROC})

set(INCLUDEDIR ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(dftbplus PROPERTIES Fortran_MODULE_DIRECTORY ${INCLUDEDIR})

target_include_directories(dftbplus INTERFACE
  $<TARGET_PROPERTY:xmlf90,INTERFACE_INCLUDE_DIRECTORIES>
  ${INCLUDEDIR})

set(OTHERLIBS xmlf90)

if (WITH_SOCKETS)
  target_include_directories(dftbplus INTERFACE
    $<TARGET_PROPERTY:fsockets,INTERFACE_INCLUDE_DIRECTORIES>)
  list(APPEND OTHERLIBS fsockets)
endif()

if(WITH_DFTD3)
  target_include_directories(dftbplus INTERFACE
    $<TARGET_PROPERTY:dftd3,INTERFACE_INCLUDE_DIRECTORIES>)
  list(APPEND OTHERLIBS dftd3)
endif()

if(WITH_MPI)
  target_include_directories(dftbplus INTERFACE
    $<TARGET_PROPERTY:mpifx,INTERFACE_INCLUDE_DIRECTORIES>)
  target_include_directories(dftbplus INTERFACE
    $<TARGET_PROPERTY:scalapackfx,INTERFACE_INCLUDE_DIRECTORIES>)
  list(APPEND OTHERLIBS mpifx scalapackfx ${SCALAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES})
endif()

target_link_libraries(dftbplus ${OTHERLIBS})

_register_install_target(dftbplus)
_register_install_modfile_dirs(${INCLUDEDIR})

get_target_property(moddirs xmlf90 INTERFACE_INCLUDE_DIRECTORIES)
_register_install_modfile_dirs(${moddirs})

if (WITH_SOCKETS)
  get_target_property(moddirs fsockets INTERFACE_INCLUDE_DIRECTORIES)
  _register_install_modfile_dirs(${moddirs})
endif()

if (WITH_DFTD3)
  get_target_property(moddirs dftd3 INTERFACE_INCLUDE_DIRECTORIES)
  _register_install_modfile_dirs(${moddirs})
endif()


#
# Compile and install program
#

set(ALL-SOURCES-F90)
set(ALL-SOURCES-FPP)

add_subdirectory(prg_dftb)

_register_custom_preprocessing(
  "${FYPP}" "${FYPP-OPTS}" "F90" "f90" "${ALL-SOURCES-FPP}" ALL-SOURCES-F90-PREPROC)

add_executable("dftb+" ${ALL-SOURCES-F90} ${ALL-SOURCES-F90-PREPROC})

target_link_libraries(dftb+ dftbplus ${LAPACK_LIBRARIES} ${ARPACK_LIBRARIES})

_register_install_target(dftb+)
