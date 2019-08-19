set(sources
  fsockets.f90
  sockets.c)

add_library(fsockets_objlib OBJECT ${sources})

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(fsockets_objlib PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})

target_include_directories(fsockets_objlib PUBLIC ${includedir})
