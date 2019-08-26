set(builddir ${CMAKE_CURRENT_BINARY_DIR})

include(ExternalProject)
ExternalProject_Add(
  DftD3Project
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/origin
  BINARY_DIR ${builddir}
  CMAKE_ARGS
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
    -DLIBRARY_ONLY=TRUE
    -DINSTALL_INCLUDE_FILES=${INSTALL_INCLUDE_FILES}
    -DINSTALL_LIB_DIR=${INSTALL_LIB_DIR}
    -DINSTALL_MOD_DIR=${INSTALL_MOD_DIR}
    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
  INSTALL_COMMAND ""
  EXCLUDE_FROM_ALL TRUE)

add_library(dftd3 INTERFACE)
add_dependencies(dftd3 DftD3Project)

if(BUILD_SHARED_LIBS)
  target_link_libraries(dftd3 INTERFACE
    $<BUILD_INTERFACE:${builddir}/lib/libdftd3.so>
    $<INSTALL_INTERFACE:${INSTALL_LIB_DIR}/libdftd3.so>)
else()
  target_link_libraries(dftd3 INTERFACE
    $<BUILD_INTERFACE:${builddir}/lib/libdftd3.a>
    $<INSTALL_INTERFACE:${INSTALL_LIB_DIR}/libdftd3.a>)
endif()

target_include_directories(dftd3 INTERFACE 
  $<BUILD_INTERFACE:${builddir}/lib/include>
  $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>)

install(TARGETS dftd3 EXPORT dftbplus-targets)


set(EXPORTED_COMPILED_LIBRARIES ${EXPORTED_COMPILED_LIBRARIES} dftd3 PARENT_SCOPE)
