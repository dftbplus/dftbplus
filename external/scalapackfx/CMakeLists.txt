set(builddir ${CMAKE_CURRENT_BINARY_DIR})

include(ExternalProject)
ExternalProject_Add(
  ScalapackfxProject
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/origin
  BINARY_DIR ${builddir}
  CMAKE_ARGS
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
    -DLIBRARY_ONLY=TRUE
    -DINSTALL_LIB_DIR=${INSTALL_LIB_DIR}
    -DINSTALL_MOD_DIR=${INSTALL_MOD_DIR}
  INSTALL_COMMAND ""
  EXCLUDE_FROM_ALL TRUE)

add_library(scalapackfx INTERFACE)
add_dependencies(scalapackfx ScalapackfxProject)

target_link_libraries(scalapackfx INTERFACE
  $<BUILD_INTERFACE:${builddir}/src/libscalapackfx.a>
  $<INSTALL_INTERFACE:${INSTALL_LIB_DIR}/libscalapackfx.a>)

target_include_directories(scalapackfx INTERFACE 
  $<BUILD_INTERFACE:${builddir}/src/include>
  $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>)

install(TARGETS scalapackfx EXPORT dftbplus-targets)

set(EXPORTED_COMPILED_LIBRARIES ${EXPORTED_COMPILED_LIBRARIES} scalapackfx PARENT_SCOPE)
