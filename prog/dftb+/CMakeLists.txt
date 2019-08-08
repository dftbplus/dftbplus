#
# General options for all targets
#

dftbp_get_fypp_defines(fypp-opts)
list(APPEND fypp-opts -I${CMAKE_CURRENT_SOURCE_DIR}/include -DRELEASE="'${RELEASE}'")


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

dftbp_register_preprocessing(
  "${FYPP}" "${fypp-opts}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)

set(dftbplus-components)
set(dftbplus-extlibs)

# Library components

list(APPEND dftbplus-components ${ALL-SOURCES-F90} ${all-sources-f90-preproc})

if(MONOLITHIC_LIBDFTBPLUS)

  list(APPEND dftbplus-components $<TARGET_OBJECTS:xmlf90_objlib>)

  if(WITH_SOCKETS)
    list(APPEND dftbplus-components $<TARGET_OBJECTS:fsockets_objlib>)
  endif()

  if(WITH_MPI)
    list(APPEND dftbplus-components $<TARGET_OBJECTS:mpifx_objlib>)
    list(APPEND dftbplus-components $<TARGET_OBJECTS:scalapackfx_objlib>)
  endif()
  
endif()

add_library(dftbplus ${dftbplus-components})

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(dftbplus PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})

# External libraries

if(NOT MONOLITHIC_LIBDFTBPLUS)
  list(APPEND dftbplus-extlibs xmlf90)
endif()

target_include_directories(dftbplus PUBLIC
  $<TARGET_PROPERTY:xmlf90,INTERFACE_INCLUDE_DIRECTORIES>
  ${includedir})

if (WITH_SOCKETS)
  if(NOT MONOLITHIC_LIBDFTBPLUS)
    list(APPEND dftbplus-extlibs fsockets)
  endif()
  target_include_directories(dftbplus PUBLIC
    $<TARGET_PROPERTY:fsockets,INTERFACE_INCLUDE_DIRECTORIES>)
endif()

if(WITH_DFTD3)
  target_include_directories(dftbplus PUBLIC
    $<TARGET_PROPERTY:dftd3,INTERFACE_INCLUDE_DIRECTORIES>)
  list(APPEND dftbplus-extlibs dftd3)
endif()

if(WITH_MPI)
  target_include_directories(dftbplus PUBLIC
    $<TARGET_PROPERTY:mpifx,INTERFACE_INCLUDE_DIRECTORIES>)
  target_include_directories(dftbplus PUBLIC
    $<TARGET_PROPERTY:scalapackfx_objlib,INTERFACE_INCLUDE_DIRECTORIES>)
  if(NOT MONOLITHIC_LIBDFTBPLUS)
    list(APPEND dftbplus-extlibs mpifx scalapackfx)
  endif()
  list(APPEND dftbplus-extlibs ${SCALAPACK_LIBRARIES})
endif()

target_link_libraries(dftbplus ${dftbplus-extlibs})

dftbp_register_install_target(dftbplus)
dftbp_register_install_mod_dirs(${includedir})

get_target_property(moddirs xmlf90 INTERFACE_INCLUDE_DIRECTORIES)
dftbp_register_install_mod_dirs(${moddirs})

if (WITH_SOCKETS)
  get_target_property(moddirs fsockets INTERFACE_INCLUDE_DIRECTORIES)
  dftbp_register_install_mod_dirs(${moddirs})
endif()

if (WITH_DFTD3)
  get_target_property(moddirs dftd3 INTERFACE_INCLUDE_DIRECTORIES)
  dftbp_register_install_mod_dirs(${moddirs})
endif()


#
# Compile and install program
#

set(ALL-SOURCES-F90)
set(ALL-SOURCES-FPP)

add_subdirectory(prg_dftb)

dftbp_register_preprocessing(
  "${FYPP}" "${fypp-opts}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)

add_executable("dftb+" ${ALL-SOURCES-F90} ${all-sources-f90-preproc})

target_link_libraries(dftb+ dftbplus ${LAPACK_LIBRARIES} ${ARPACK_LIBRARIES})

dftbp_register_install_target(dftb+)
