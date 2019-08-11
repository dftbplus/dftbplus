#
# General options for all targets
#
list(APPEND FYPP_FLAGS -I${CMAKE_CURRENT_SOURCE_DIR}/include -DRELEASE="'${RELEASE}'")


#
# Compile and install library
#

set(ALL-SOURCES-F90)
set(ALL-SOURCES-FPP)

add_subdirectory(lib_common)
add_subdirectory(lib_derivs)
add_subdirectory(lib_dftb)
add_subdirectory(lib_dftbplus)
add_subdirectory(lib_elecsolvers)
add_subdirectory(lib_extlibs)
add_subdirectory(lib_geoopt)
add_subdirectory(lib_io)
add_subdirectory(lib_math)
add_subdirectory(lib_md)
add_subdirectory(lib_mixer)
add_subdirectory(lib_timedep)
add_subdirectory(lib_type)
if(BUILD_API)
  add_subdirectory(api/mm)
endif()

dftbp_register_preprocessing(
  "${FYPP}" "${FYPP_FLAGS}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)

set(dftbplus-components)
#set(dftbplus-extlibs)
set(dftbplus-includedirs)

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
target_include_directories(dftbplus PUBLIC ${includedir})
if(BUILD_API)
  target_include_directories(dftbplus INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/api/mm)
endif()

# External libraries

if(MONOLITHIC_LIBDFTBPLUS)
  target_include_directories(dftbplus PUBLIC
    $<TARGET_PROPERTY:xmlf90,INTERFACE_INCLUDE_DIRECTORIES>
    ${includedir})
else()
  target_link_libraries(dftbplus PUBLIC xmlf90)
endif()

if(WITH_SOCKETS)
  if(MONOLITHIC_LIBDFTBPLUS)
    target_include_directories(dftbplus PUBLIC
      $<TARGET_PROPERTY:fsockets,INTERFACE_INCLUDE_DIRECTORIES>)
  else()
    target_link_libraries(dftbplus PUBLIC fsockets)
  endif()
endif()

if(WITH_DFTD3)
  target_link_libraries(dftbplus PUBLIC dftd3)
endif()

if(WITH_TRANSPORT)
  target_link_libraries(dftbplus PUBLIC libnegf poisson)
endif()

if(WITH_ELSI)
  target_link_libraries(dftbplus PUBLIC ${ELSI_LIBRARIES})
  target_include_directories(dftbplus PUBLIC ${ELSI_INCLUDE_DIRS})
endif()

if(WITH_PEXSI)
  target_link_libraries(dftbplus PUBLIC ${PEXSI_LIBRARIES})
  target_include_directories(dftbplus PUBLIC ${PEXSI_INCLUDE_DIRS})
endif()

if(WITH_MPI)
  if(MONOLITHIC_LIBDFTBPLUS)
    target_include_directories(dftbplus PUBLIC
      $<TARGET_PROPERTY:mpifx,INTERFACE_INCLUDE_DIRECTORIES>)
    target_include_directories(dftbplus PUBLIC
      $<TARGET_PROPERTY:scalapackfx_objlib,INTERFACE_INCLUDE_DIRECTORIES>)
  else()
    target_link_libraries(dftbplus PUBLIC mpifx scalapackfx)
  endif()
  target_link_libraries(dftbplus PUBLIC ${SCALAPACK_LIBRARIES})
endif()

target_link_libraries(dftbplus PUBLIC ${LAPACK_LIBRARIES})

if(WITH_ARPACK)
  target_link_libraries(dftbplus PUBLIC ${ARPACK_LIBRARIES})
endif()

target_link_libraries(dftbplus PUBLIC ${OTHER_LIBRARIES})

target_include_directories(dftbplus PUBLIC ${OTHER_INCLUDE_DIRS})

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
  "${FYPP}" "${FYPP_FLAGS}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)

add_executable("dftb+" ${ALL-SOURCES-F90} ${all-sources-f90-preproc})

target_link_libraries(dftb+ dftbplus)

dftbp_register_install_target(dftb+)
