set(builddir ${CMAKE_CURRENT_BINARY_DIR})

include(ExternalProject)
ExternalProject_Add(
  MpifxProject
  DOWNLOAD_COMMAND ""
  UPDATE_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/origin
  BINARY_DIR ${builddir}
  CMAKE_ARGS
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
    -DLIBRARY_ONLY=TRUE
    -DINSTALL_LIB_DIR=${INSTALL_LIB_DIR}
    -DINSTALL_MOD_DIR=${INSTALL_MOD_DIR}
  INSTALL_COMMAND cmake -E echo
  EXCLUDE_FROM_ALL TRUE)

add_library(mpifx INTERFACE)
add_dependencies(mpifx MpifxProject)

target_link_libraries(mpifx INTERFACE
  $<BUILD_INTERFACE:${builddir}/src/libmpifx.a>
  $<INSTALL_INTERFACE:${INSTALL_LIB_DIR}/libmpifx.a>)

target_include_directories(mpifx INTERFACE 
  $<BUILD_INTERFACE:${builddir}/src/include>
  $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>)

install(TARGETS mpifx EXPORT dftbplus-targets)
