# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

add_subdirectory(${TBLITE_SOURCE_DIR} ${TBLITE_BINARY_DIR})

add_library(tblite::tblite INTERFACE IMPORTED GLOBAL)
target_link_libraries(tblite::tblite INTERFACE tblite)
