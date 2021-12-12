# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

add_subdirectory(${TOML-F_SOURCE_DIR} ${TOML-F_BINARY_DIR})

add_library(toml-f::toml-f INTERFACE IMPORTED GLOBAL)
target_link_libraries(toml-f::toml-f INTERFACE toml-f)
