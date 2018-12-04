set(SOURCES
  fsockets.f90
  sockets.c)

add_library(fsockets ${SOURCES})

set(INCLUDEDIR ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(fsockets PROPERTIES Fortran_MODULE_DIRECTORY ${INCLUDEDIR})

target_include_directories(fsockets INTERFACE ${INCLUDEDIR})
