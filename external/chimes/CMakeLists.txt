# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!
# Note: ChIMES currently can only be invoked as a submodule!

# Directly build ChIMES as its CMake builder is currently rather inclusion unfriendly
set(sources
  origin/chimesFF/src/chimesFF.cpp
  origin/serial_interface/src/serial_chimes_interface.cpp
  origin/serial_interface/api/wrapper-C.cpp
  origin/serial_interface/api/wrapper-F.F90
  origin/serial_interface/api/wrapper-F08.F90
)

add_library(chimes ${sources})
target_include_directories(chimes
  PRIVATE origin/chimesFF/src origin/serial_interface/src origin/serial_interface/api)

set(includedir ${CMAKE_CURRENT_BINARY_DIR}/include)
set_target_properties(chimes PROPERTIES Fortran_MODULE_DIRECTORY ${includedir})
set_target_properties(chimes PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})
target_include_directories(chimes PUBLIC ${includedir})

add_library(Chimes::Chimes INTERFACE IMPORTED GLOBAL)
target_link_libraries(Chimes::Chimes INTERFACE chimes)
