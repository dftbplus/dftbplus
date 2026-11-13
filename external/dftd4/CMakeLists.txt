# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Propagate OpenMP option correctly to subproject
set(DFTD4_WITH_OpenMP ${WITH_OMP})

# Do not compile with C-API support as we do not want to vendor the subproject C-API
set(DFTD4_WITH_API FALSE)
# Subproject uses both cache variables to steer API build
set(WITH_API FALSE)

add_subdirectory(${DFTD4_SOURCE_DIR} ${DFTD4_BINARY_DIR})
