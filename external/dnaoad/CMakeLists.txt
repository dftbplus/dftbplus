# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

set(DNAOAD_INSTALL ON)
set(DNAOAD_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
add_subdirectory(${DNAOAD_SOURCE_DIR} ${DNAOAD_BINARY_DIR})
