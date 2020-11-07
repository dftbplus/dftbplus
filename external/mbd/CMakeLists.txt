# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

set(ENABLE_SCALAPACK_MPI ${WITH_MPI})
set(ENABLE_ELSI ${WITH_ELSI})
set(ENABLE_C_API OFF)
set(SCALAPACK_LIBRARIES "${SCALAPACK_LIBRARY}")
set(MKL_LIBRARIES "")
set(BUILD_TESTING FALSE)

add_subdirectory(${MBD_SOURCE_DIR} ${MBD_BINARY_DIR})

add_library(Mbd::Mbd INTERFACE IMPORTED GLOBAL)
target_link_libraries(Mbd::Mbd INTERFACE Mbd)
