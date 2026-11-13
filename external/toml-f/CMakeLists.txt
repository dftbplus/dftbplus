# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Disable compilation of unit tests (would also trigger fetching test-drive)
# (subproject requires an old CMake version, so option() must be told to honor normal variables)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(BUILD_TESTING FALSE)

add_subdirectory(${TOML-F_SOURCE_DIR} ${TOML-F_BINARY_DIR})
