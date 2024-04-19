# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

set(BUILD_TESTING FALSE)

get_target_property(LAPACK_LIBRARIES LAPACK::LAPACK INTERFACE_LINK_LIBRARIES)

# This setting is required by libNEGF but never really used.
set(LAPACK_LIBRARY_DIRS "dummy")

add_subdirectory(${NEGF_SOURCE_DIR} ${NEGF_BINARY_DIR})

add_library(Negf::Negf INTERFACE IMPORTED GLOBAL)
target_link_libraries(Negf::Negf INTERFACE Negf)
