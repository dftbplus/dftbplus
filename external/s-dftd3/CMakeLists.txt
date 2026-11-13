# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Propagate OpenMP option correctly to subproject
set(SDFTD3_WITH_OpenMP ${WITH_OMP})

# Do not compile with C-API support as we do not want to vendor the subproject C-API
set(SDFTD3_WITH_API FALSE)

# Disable MPI
set(SDFTD3_WITH_MPI FALSE)

# Disable compilation of unit tests
set(SDFTD3_WITH_TESTS FALSE)

add_subdirectory(${S-DFTD3_SOURCE_DIR} ${S-DFTD3_BINARY_DIR})
