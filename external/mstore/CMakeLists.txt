# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Propagate OpenMP option correctly to subproject
set(WITH_OpenMP ${WITH_OMP})

add_subdirectory(${MSTORE_SOURCE_DIR} ${MSTORE_BINARY_DIR})

add_library(mstore::mstore INTERFACE IMPORTED GLOBAL)
target_link_libraries(mstore::mstore INTERFACE mstore)
