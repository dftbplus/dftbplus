#
# General options for all targets
#
list(APPEND FYPP_FLAGS -I${CMAKE_CURRENT_SOURCE_DIR}/include -DRELEASE="'${RELEASE}'")

#
# Compile and install library
#

set(ALL-SOURCES-F90)
set(ALL-SOURCES-FPP)
set(ALL-SOURCES-C)

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
if(WITH_API)
  add_subdirectory(api/mm)
endif()
if(WITH_TRANSPORT)
  add_subdirectory(lib_poisson)
endif()
if(WITH_GPU)
  add_subdirectory(lib_magmahelper)
endif()

dftbp_preprocess("${FYPP}" "${FYPP_FLAGS}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)


# Library components

add_library(dftbplus ${ALL-SOURCES-F90} ${all-sources-f90-preproc} ${ALL-SOURCES-C})

target_sources(dftbplus PRIVATE $<TARGET_OBJECTS:xmlf90_objlib>)

if(WITH_SOCKETS)
  target_sources(dftbplus PRIVATE $<TARGET_OBJECTS:fsockets_objlib>)
endif()

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(dftbplus PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})
target_include_directories(dftbplus PUBLIC
  $<BUILD_INTERFACE:${includedir}> $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>)
if(WITH_API)
  target_include_directories(dftbplus INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/api/mm> $<INSTALL_INTERFACE:${INSTALL_INC_DIR}>)
endif()


# External libraries

target_include_directories(dftbplus PUBLIC
  $<BUILD_INTERFACE:$<TARGET_PROPERTY:xmlf90_objlib,INTERFACE_INCLUDE_DIRECTORIES>>)

target_include_directories(dftbplus PRIVATE ${CMAKE_SOURCE_DIR}/external/dftd4refs)

if(WITH_SOCKETS)
  target_include_directories(dftbplus PUBLIC
    $<BUILD_INTERFACE:$<TARGET_PROPERTY:fsockets_objlib,INTERFACE_INCLUDE_DIRECTORIES>>)
endif()

if(WITH_DFTD3)
  target_link_libraries(dftbplus PUBLIC dftd3)
endif()

if(WITH_TRANSPORT)
  target_link_libraries(dftbplus PUBLIC negf sparskit mudpack)
endif()

if(WITH_ELSI)
  # Some compilers need all mod files (also indirect ones), so include dir must be exported
  target_link_libraries(dftbplus PUBLIC elsi::elsi)
  target_link_libraries(dftbplus PRIVATE ${ELSI_EXTERNAL_LIBRARIES})
  if(ELSI_WITH_PEXSI)
    target_link_libraries(dftbplus PRIVATE ${PEXSI_EXTERNAL_LIBRARIES})
  endif()
endif()

if(WITH_PLUMED)
  #target_link_libraries(dftbplus PRIVATE PkgConfig::PLUMED)
  target_link_libraries(dftbplus PRIVATE ${PLUMED_LIBRARIES})
endif()

if(WITH_MPI)
  target_link_libraries(dftbplus PUBLIC mpifx scalapackfx)
  target_link_libraries(dftbplus PRIVATE ${SCALAPACK_LIBRARIES})
endif()

target_link_libraries(dftbplus PRIVATE ${LAPACK_LIBRARIES})

if(WITH_ARPACK)
  target_link_libraries(dftbplus PRIVATE ${ARPACK_LIBRARIES})
endif()

if (WITH_GPU)
  target_include_directories(dftbplus PUBLIC ${MAGMA_INCLUDE_DIRS})
  target_link_libraries(dftbplus PRIVATE ${MAGMA_LIBRARIES})
endif()

target_link_libraries(dftbplus PRIVATE ${OTHER_LIBRARIES})

target_include_directories(dftbplus PUBLIC ${OTHER_INCLUDE_DIRS})


# Installation

if(WITH_API OR BUILD_SHARED_LIBS)
  install(TARGETS dftbplus DESTINATION ${INSTALL_LIB_DIR} EXPORT dftbplus-targets)
endif()

if(INSTALL_INCLUDE_FILES)
  install(DIRECTORY ${includedir}/ DESTINATION ${INSTALL_MOD_DIR})
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/api/mm/dftbplus.h DESTINATION ${INSTALL_INC_DIR})

  get_target_property(moddirs xmlf90_objlib INTERFACE_INCLUDE_DIRECTORIES)
  foreach(moddir IN LISTS moddirs)
    install(DIRECTORY "${moddir}/" DESTINATION "${INSTALL_MOD_DIR}")
  endforeach()

  if (WITH_SOCKETS)
    get_target_property(moddirs fsockets_objlib INTERFACE_INCLUDE_DIRECTORIES)
    foreach(moddir IN LISTS moddirs)
      install(DIRECTORY "${moddir}/" DESTINATION "${INSTALL_MOD_DIR}")
    endforeach()
  endif()
endif()


#
# Compile and install program
#

set(ALL-SOURCES-F90)
set(ALL-SOURCES-FPP)

add_subdirectory(prg_dftb)

dftbp_preprocess("${FYPP}" "${FYPP_FLAGS}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)

add_executable(dftb+ ${ALL-SOURCES-F90} ${all-sources-f90-preproc})

target_link_libraries(dftb+ dftbplus)

install(TARGETS dftb+ DESTINATION ${INSTALL_BIN_DIR})
