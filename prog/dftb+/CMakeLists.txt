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

dftbp_preprocess("${FYPP}" "${FYPP_FLAGS}" "F90" "f90" "${ALL-SOURCES-FPP}" all-sources-f90-preproc)


# Library components

add_library(dftbplus ${ALL-SOURCES-F90} ${all-sources-f90-preproc})

target_sources(dftbplus PRIVATE $<TARGET_OBJECTS:xmlf90_objlib>)

if(WITH_SOCKETS)
  target_sources(dftbplus PRIVATE $<TARGET_OBJECTS:fsockets_objlib>)
endif()

if(WITH_GPU)
  target_sources(dftbplus PRIVATE $<TARGET_OBJECTS:magmahelper_objlib>)
endif()

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(dftbplus PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})
target_include_directories(dftbplus PUBLIC ${includedir})
if(BUILD_API)
  target_include_directories(dftbplus INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/api/mm)
endif()


# External libraries

target_include_directories(dftbplus PUBLIC
  $<TARGET_PROPERTY:xmlf90_objlib,INTERFACE_INCLUDE_DIRECTORIES> ${includedir})

if(WITH_SOCKETS)
  target_include_directories(dftbplus PUBLIC
    $<TARGET_PROPERTY:fsockets_objlib,INTERFACE_INCLUDE_DIRECTORIES>)
endif()

if(WITH_GPU)
  target_include_directories(dftbplus PUBLIC
    $<TARGET_PROPERTY:magmahelper_objlib,INTERFACE_INCLUDE_DIRECTORIES>)
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
  target_link_libraries(dftbplus PUBLIC mpifx scalapackfx)
  target_link_libraries(dftbplus PUBLIC ${SCALAPACK_LIBRARIES})
endif()

target_link_libraries(dftbplus PUBLIC ${LAPACK_LIBRARIES})

if(WITH_ARPACK)
  target_link_libraries(dftbplus PUBLIC ${ARPACK_LIBRARIES})
endif()

target_link_libraries(dftbplus PUBLIC ${OTHER_LIBRARIES})

target_include_directories(dftbplus PUBLIC ${OTHER_INCLUDE_DIRS})


# Installation

if(BUILD_API)
  install(TARGETS dftbplus DESTINATION ${INSTALL_LIB_DIR})
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
  
  if (WITH_GPU)
    get_target_property(moddirs magmahelper_objlib INTERFACE_INCLUDE_DIRECTORIES)
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

add_executable("dftb+" ${ALL-SOURCES-F90} ${all-sources-f90-preproc})

target_link_libraries(dftb+ dftbplus)

dftbp_install_target(dftb+)
