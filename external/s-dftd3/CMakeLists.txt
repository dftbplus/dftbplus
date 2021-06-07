# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

add_subdirectory(${S-DFTD3_SOURCE_DIR} ${S-DFTD3_BINARY_DIR})

add_library(s-dftd3::s-dftd3 INTERFACE IMPORTED GLOBAL)
target_link_libraries(s-dftd3::s-dftd3 INTERFACE s-dftd3)
