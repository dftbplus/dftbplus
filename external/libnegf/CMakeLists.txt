# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

set(BUILD_TESTING FALSE)

add_subdirectory(${NEGF_SOURCE_DIR} ${NEGF_BINARY_DIR})

add_library(Negf::Negf INTERFACE IMPORTED GLOBAL)
target_link_libraries(Negf::Negf INTERFACE Negf)
