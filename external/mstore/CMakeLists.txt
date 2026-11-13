# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

add_subdirectory(${MSTORE_SOURCE_DIR} ${MSTORE_BINARY_DIR})

add_library(mstore::mstore INTERFACE IMPORTED GLOBAL)
target_link_libraries(mstore::mstore INTERFACE mstore)
