# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Propagate OpenMP option correctly to subproject
set(WITH_OpenMP ${WITH_OMP})

# Do not compile with C-API support as we do not want to vendor the subproject C-API
set(WITH_API FALSE)

# We have an optimized BLAS implementation here, therefore it should be used
set(WITH_BLAS TRUE)

add_subdirectory(${S-DFTD3_SOURCE_DIR} ${S-DFTD3_BINARY_DIR})

add_library(s-dftd3::s-dftd3 INTERFACE IMPORTED GLOBAL)
target_link_libraries(s-dftd3::s-dftd3 INTERFACE s-dftd3)
