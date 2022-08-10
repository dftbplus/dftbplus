set(projectdir ${PROJECT_SOURCE_DIR})

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

add_subdirectory(common)
add_subdirectory(derivs)
add_subdirectory(dftb)
add_subdirectory(dftbplus)
add_subdirectory(elecsolvers)
add_subdirectory(extlibs)
add_subdirectory(geoopt)
add_subdirectory(io)
add_subdirectory(math)
add_subdirectory(md)
add_subdirectory(mixer)
add_subdirectory(modelindependent)
add_subdirectory(reks)
add_subdirectory(timedep)
add_subdirectory(type)
add_subdirectory(solvation)
add_subdirectory(transport)
if(WITH_API)
  add_subdirectory(api/mm)
endif()
if(WITH_POISSON)
  add_subdirectory(poisson)
endif()

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)

# Library components

add_library(dftbplus ${ALL-SOURCES-F90} ${all-sources-f90-preproc} ${ALL-SOURCES-C})

target_sources(dftbplus PRIVATE $<TARGET_OBJECTS:xmlf90_objlib>)

# Link ddCOSMO implementation
target_sources(dftbplus PRIVATE $<TARGET_OBJECTS:ddcosmo_objlib>)

if(WITH_SOCKETS)
  target_sources(dftbplus PRIVATE $<TARGET_OBJECTS:fsockets_objlib>)
endif()

if(EXTERNAL_MODEL_INTERNAL)
  target_sources(dftbplus PRIVATE $<TARGET_OBJECTS:externalmodel_objlib>)
else()

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

target_include_directories(dftbplus PUBLIC
  $<BUILD_INTERFACE:$<TARGET_PROPERTY:ddcosmo_objlib,INTERFACE_INCLUDE_DIRECTORIES>>)

target_include_directories(dftbplus PRIVATE ${projectdir}/external/dftd4refs)

if(WITH_SOCKETS)
  target_include_directories(dftbplus PUBLIC
    $<BUILD_INTERFACE:$<TARGET_PROPERTY:fsockets_objlib,INTERFACE_INCLUDE_DIRECTORIES>>)
endif()

if(EXTERNAL_MODEL_INTERNAL)
  target_include_directories(dftbplus PUBLIC
    $<BUILD_INTERFACE:$<TARGET_PROPERTY:externalmodel_objlib,INTERFACE_INCLUDE_DIRECTORIES>>)
endif()

if(WITH_MBD)
  target_link_libraries(dftbplus PUBLIC Mbd::Mbd)
endif()

if(WITH_TRANSPORT)
  target_link_libraries(dftbplus PUBLIC Negf::Negf)
endif()

if(WITH_POISSON)
  target_link_libraries(dftbplus PUBLIC mudpack)
endif()

if(WITH_SDFTD3)
  target_link_libraries(dftbplus PUBLIC s-dftd3::s-dftd3)
endif()
if(WITH_TBLITE)
  target_link_libraries(dftbplus PUBLIC tblite::tblite)
endif()
if(WITH_TBLITE OR WITH_SDFTD3)
  target_link_libraries(dftbplus PUBLIC mctc-lib::mctc-lib)
endif()

if(WITH_ELSI)
  # Some compilers need all mod files (also indirect ones), so include dir must be exported
  target_link_libraries(dftbplus PUBLIC elsi::elsi)
  target_link_libraries(dftbplus PRIVATE ${ELSI_EXTERNAL_LIBRARIES})
  if(ELSI_WITH_PEXSI)
    target_link_libraries(dftbplus PRIVATE ${PEXSI_EXTERNAL_LIBRARIES})
  endif()
  if(WITH_GPU)
    target_link_libraries(dftbplus PRIVATE CUDA::cudart)
    target_link_libraries(dftbplus PRIVATE CUDA::cublas)
  endif()
endif()

if(WITH_PLUMED)
  target_link_libraries(dftbplus PRIVATE Plumed::Plumed)
endif()

if(WITH_CHIMES)
  target_link_libraries(dftbplus PUBLIC ChimesCalc::ChimesCalc_Fortran)
endif()

if(WITH_MPI)
  target_link_libraries(dftbplus PRIVATE MPI::MPI_Fortran)
  target_link_libraries(dftbplus PUBLIC MpiFx::MpiFx ScalapackFx::ScalapackFx)
endif()

# OpenMP dependency must be public in order to forward eventual linker flags
# (workaround for NAG compiler) to targets linking to it.
target_link_libraries(dftbplus PUBLIC OpenMP::OpenMP_Fortran)

if(WITH_ARPACK)
  target_link_libraries(dftbplus PUBLIC ARPACK::ARPACK)
  if(WITH_MPI)
    target_link_libraries(dftbplus PUBLIC PARPACK::PARPACK)
  endif()
endif()

target_link_libraries(dftbplus PUBLIC LAPACK::LAPACK)

if(WITH_MAGMA)
  target_link_libraries(dftbplus PUBLIC Magma::Magma)
endif()


# Installation

if(WITH_API OR BUILD_SHARED_LIBS)
  install(TARGETS dftbplus DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    EXPORT dftbplus-targets)
endif()

if(INSTALL_INCLUDE_FILES)
  install(DIRECTORY ${includedir}/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODULEDIR}")

  if(WITH_API)
    CONFIGURE_FILE(
      "${CMAKE_CURRENT_SOURCE_DIR}/api/mm/dftbplus.h"
      "${CMAKE_CURRENT_BINARY_DIR}/dftbplus.h"
      @ONLY
      )
    install(
      FILES "${CMAKE_CURRENT_BINARY_DIR}/dftbplus.h"
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
      )
  endif()

  get_target_property(moddirs xmlf90_objlib INTERFACE_INCLUDE_DIRECTORIES)
  foreach(moddir IN LISTS moddirs)
    install(DIRECTORY "${moddir}/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODULEDIR}")
  endforeach()

  get_target_property(moddirs ddcosmo_objlib INTERFACE_INCLUDE_DIRECTORIES)
  foreach(moddir IN LISTS moddirs)
    install(DIRECTORY "${moddir}/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODULEDIR}")
  endforeach()


  if (WITH_SOCKETS)
    get_target_property(moddirs fsockets_objlib INTERFACE_INCLUDE_DIRECTORIES)
    foreach(moddir IN LISTS moddirs)
      install(DIRECTORY "${moddir}/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_MODULEDIR}")
    endforeach()
  endif()

  if (EXTERNAL_MODEL_INTERNAL)
    get_target_property(moddirs externalmodel_objlib INTERFACE_INCLUDE_DIRECTORIES)
  endif()

endif()

list(APPEND PKG_CONFIG_LIBS dftbplus)
set(PKG_CONFIG_LIBS "${PKG_CONFIG_LIBS}" PARENT_SCOPE)
