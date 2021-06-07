# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

add_subdirectory(${DFTD4_SOURCE_DIR} ${DFTD4_BINARY_DIR})

add_library(dftd4::dftd4 INTERFACE IMPORTED GLOBAL)
target_link_libraries(dftd4::dftd4 INTERFACE dftd4)
