# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

set(ENABLE_SCALAPACK_MPI ${WITH_MPI})
set(ENABLE_ELSI ${WITH_ELSI})
set(ENABLE_C_API OFF)
if(WITH_MPI)
  get_target_property(SCALAPACK_LIBRARIES Scalapack::Scalapack INTERFACE_LINK_LIBRARIES)
endif()
if(WITH_ELSI)
  get_target_property(ELSI_LIBRARIES elsi::elsi INTERFACE_LINK_LIBRARIES)
  # FIXME: Hack needed until https://github.com/libmbd/libmbd/pull/30 is resolved
  get_target_property(ELSI_INCLUDE_DIRECTORIES elsi::elsi INTERFACE_INCLUDE_DIRECTORIES)
  include_directories("${ELSI_INCLUDE_DIRECTORIES}")
endif()
get_target_property(LAPACK_LIBRARIES LAPACK::LAPACK INTERFACE_LINK_LIBRARIES)
set(MKL_LIBRARIES "")
set(BUILD_TESTING FALSE)

add_subdirectory(${MBD_SOURCE_DIR} ${MBD_BINARY_DIR})

add_library(Mbd::Mbd INTERFACE IMPORTED GLOBAL)
target_link_libraries(Mbd::Mbd INTERFACE Mbd)
