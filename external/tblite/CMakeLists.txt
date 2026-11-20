# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Propagate paralellization options
set(TBLITE_WITH_OpenMP ${WITH_OMP})
set(TBLITE_WITH_MPI FALSE)

# Disable various option which are not needed when building TBLite as a subproject
set(TBLITE_WITH_TESTS FALSE)
set(TBLITE_WITH_API FALSE)
set(TBLITE_WITH_PYTHON FALSE)

add_subdirectory(${TBLITE_SOURCE_DIR} ${TBLITE_BINARY_DIR})

target_compile_options(
  "tblite-exe"
  PRIVATE
  "$<$<Fortran_COMPILER_ID:GNU>:-fall-intrinsics>"
)
