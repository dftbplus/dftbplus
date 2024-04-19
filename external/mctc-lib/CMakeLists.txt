# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

# Propagate OpenMP option correctly to subproject
set(WITH_OpenMP ${WITH_OMP})
# Disable JSON Fortran dependency in case of subproject include
set(WITH_JSON FALSE)

add_subdirectory(${MCTC-LIB_SOURCE_DIR} ${MCTC-LIB_BINARY_DIR})

add_library(mctc-lib::mctc-lib INTERFACE IMPORTED GLOBAL)
target_link_libraries(mctc-lib::mctc-lib INTERFACE mctc-lib)
