# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Propagate OpenMP option correctly to subproject
set(WITH_OpenMP ${WITH_OMP})

# Disable compilation of unit tests
set(WITH_TESTS FALSE)

# Do not compile with C-API support as we do not want to vendor the subproject C-API
set(WITH_API FALSE)

add_subdirectory(${TBLITE_SOURCE_DIR} ${TBLITE_BINARY_DIR})

target_compile_options(
  "tblite-exe"
  PRIVATE
  "$<$<Fortran_COMPILER_ID:GNU>:-fall-intrinsics>"
)

add_library(tblite::tblite INTERFACE IMPORTED GLOBAL)
target_link_libraries(tblite::tblite INTERFACE tblite)
