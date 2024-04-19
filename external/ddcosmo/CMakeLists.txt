set(
  sources
  "ddcosmo_core.F90"
  "ddcosmo_solver.F90"
)

add_library(
  "ddcosmo_objlib"
  OBJECT
  ${sources}
)

set_target_properties(
  "ddcosmo_objlib"
  PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
set_target_properties(
  "ddcosmo_objlib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE "${BUILD_SHARED_LIBS}"
)
target_include_directories(
  "ddcosmo_objlib"
  PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/include"
)
if(WITH_OMP)
  target_link_libraries(
    "ddcosmo_objlib"
    PUBLIC OpenMP::OpenMP_Fortran
  )
endif()
