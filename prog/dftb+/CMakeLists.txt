file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../" projectdir)
#
# General options for all targets
#
set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${CMAKE_CURRENT_SOURCE_DIR}/include -DRELEASE="'${RELEASE}'")

if(WITH_API)
  message(STATUS "DFTB+ API version: ${API_VERSION}")
  list(APPEND fypp_flags -DAPIVERSION="'${API_VERSION}'")
  list(APPEND fypp_flags -DAPIMAJOR=${API_VERSION_MAJOR} -DAPIMINOR=${API_VERSION_MINOR}
    -DAPIPATCH=${API_VERSION_PATCH})
endif()

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
add_subdirectory(lib_poisson)
add_subdirectory(lib_reks)
add_subdirectory(lib_timedep)
add_subdirectory(lib_type)
add_subdirectory(lib_solvation)
if(WITH_API)
  add_subdirectory(api/mm)
endif()

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)


# Library components

add_library(dftbplus ${ALL-SOURCES-F90} ${all-sources-f90-preproc} ${ALL-SOURCES-C})

target_sources(dftbplus PRIVATE $<TARGET_OBJECTS:xmlf90_objlib>)

if(WITH_SOCKETS)
  target_sources(dftbplus PRIVATE $<TARGET_OBJECTS:fsockets_objlib>)
endif()

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(dftbplus PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})
target_include_directories(dftbplus PUBLIC
  $<BUILD_INTERFACE:${includedir}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODULEDIR}>)

if(WITH_API)
  target_include_directories(dftbplus INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/api/mm>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_INCLUDEDIR}>)
endif()


# External libraries

target_include_directories(dftbplus PUBLIC
  $<BUILD_INTERFACE:$<TARGET_PROPERTY:xmlf90_objlib,INTERFACE_INCLUDE_DIRECTORIES>>)

target_include_directories(dftbplus PRIVATE ${projectdir}/external/dftd4refs)

if(WITH_SOCKETS)
  target_include_directories(dftbplus PUBLIC
    $<BUILD_INTERFACE:$<TARGET_PROPERTY:fsockets_objlib,INTERFACE_INCLUDE_DIRECTORIES>>)
endif()

if(WITH_DFTD3)
  target_link_libraries(dftbplus PUBLIC dftd3)
endif()

if(WITH_MBD)
  target_link_libraries(dftbplus PUBLIC Mbd::Mbd)
endif()

if(WITH_TRANSPORT)
  target_link_libraries(dftbplus PUBLIC Negf::Negf)
endif()
target_link_libraries(dftbplus PUBLIC mudpack)

if(WITH_ELSI)
  # Some compilers need all mod files (also indirect ones), so include dir must be exported
  target_link_libraries(dftbplus PUBLIC elsi::elsi)
  target_link_libraries(dftbplus PRIVATE ${ELSI_EXTERNAL_LIBRARIES})
  if(ELSI_WITH_PEXSI)
    target_link_libraries(dftbplus PRIVATE ${PEXSI_EXTERNAL_LIBRARIES})
  endif()
endif()

if(WITH_PLUMED)
  target_link_libraries(dftbplus PRIVATE Plumed::Plumed)
endif()

if(WITH_MPI)
  target_link_libraries(dftbplus PUBLIC MPI::MPI_Fortran)
  target_link_libraries(dftbplus PUBLIC MpiFx::MpiFx ScalapackFx::ScalapackFx)
endif()

if(WITH_OMP)
  # OpenMP dependency must be public in order to forward eventual linker flags
  # (workaround for NAG compiler) to targets linking to it.
  target_link_libraries(dftbplus PUBLIC OpenMP::OpenMP_Fortran)
endif()

if(WITH_ARPACK)
  target_link_libraries(dftbplus PRIVATE Arpack::Arpack)
endif()

target_link_libraries(dftbplus PRIVATE LAPACK::LAPACK)

if (WITH_GPU)
  target_link_libraries(dftbplus PUBLIC Magma::Magma)
endif()

target_link_libraries(dftbplus PRIVATE ${OTHER_LIBRARIES})

target_include_directories(dftbplus PUBLIC ${OTHER_INCLUDE_DIRS})


# Installation

if(WITH_API OR BUILD_SHARED_LIBS)
  install(TARGETS dftbplus DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    EXPORT dftbplus-targets)
endif()

if(INSTALL_INCLUDE_FILES)
  install(DIRECTORY ${includedir}/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODULEDIR}")
  CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/api/mm/dftbplus.h"
    "${CMAKE_CURRENT_BINARY_DIR}/dftbplus.h"
    @ONLY
  )
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/dftbplus.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )

  get_target_property(moddirs xmlf90_objlib INTERFACE_INCLUDE_DIRECTORIES)
  foreach(moddir IN LISTS moddirs)
    install(DIRECTORY "${moddir}/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODULEDIR}")
  endforeach()

  if (WITH_SOCKETS)
    get_target_property(moddirs fsockets_objlib INTERFACE_INCLUDE_DIRECTORIES)
    foreach(moddir IN LISTS moddirs)
      install(DIRECTORY "${moddir}/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODULEDIR}")
    endforeach()
  endif()
endif()

list(APPEND PKG_CONFIG_LIBS dftbplus)
set(PKG_CONFIG_LIBS "${PKG_CONFIG_LIBS}" PARENT_SCOPE)

#
# Compile and install program
#
if(NOT BUILD_EXPORTED_TARGETS_ONLY)

  set(ALL-SOURCES-F90)
  set(ALL-SOURCES-FPP)

  add_subdirectory(prg_dftb)

  dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${ALL-SOURCES-FPP}"
    all-sources-f90-preproc)

  add_executable(dftb+ ${ALL-SOURCES-F90} ${all-sources-f90-preproc})

  target_link_libraries(dftb+ dftbplus)

  install(TARGETS dftb+ DESTINATION ${CMAKE_INSTALL_BINDIR})

endif()
