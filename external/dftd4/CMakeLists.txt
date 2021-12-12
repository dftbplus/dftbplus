# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Propagate OpenMP option correctly to subproject
set(WITH_OpenMP ${WITH_OMP})

# Do not compile with C-API support as we do not want to vendor the subproject C-API
set(WITH_API FALSE)

add_subdirectory(${DFTD4_SOURCE_DIR} ${DFTD4_BINARY_DIR})

add_library(dftd4::dftd4 INTERFACE IMPORTED GLOBAL)
target_link_libraries(dftd4::dftd4 INTERFACE dftd4)
