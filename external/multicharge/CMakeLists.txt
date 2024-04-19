# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Propagate OpenMP option correctly to subproject
set(WITH_OpenMP ${WITH_OMP})

add_subdirectory(${MULTICHARGE_SOURCE_DIR} ${MULTICHARGE_BINARY_DIR})

add_library(multicharge::multicharge INTERFACE IMPORTED GLOBAL)
target_link_libraries(multicharge::multicharge INTERFACE multicharge)
