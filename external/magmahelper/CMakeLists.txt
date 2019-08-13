set(sources
  count_gpus.c
  magma_device_info.f90)

add_library(magmahelper_objlib OBJECT ${sources})
add_library(magmahelper $<TARGET_OBJECTS:magmahelper_objlib>)

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)

set_target_properties(magmahelper_objlib PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})

target_include_directories(magmahelper_objlib PUBLIC ${includedir})
target_include_directories(magmahelper PUBLIC ${includedir})
