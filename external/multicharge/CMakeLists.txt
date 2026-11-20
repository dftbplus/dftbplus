# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Propagate OpenMP option correctly to subproject
set(MULTICHARGE_WITH_OpenMP ${WITH_OMP})

# Disable compilation of unit tests
set(MULTICHARGE_WITH_TESTS FALSE)

add_subdirectory(${MULTICHARGE_SOURCE_DIR} ${MULTICHARGE_BINARY_DIR})
