set(INSTALL_EXPORT_NAME dftbplus-targets)
set(MPIFX_INCLUDE_DIR "${CMAKE_BINARY_DIR}/external/mpifx/origin/lib/include" CACHE STRING
  "Where to find the MPIFX include directory")
# prevent libNEGF testers being built if inside DFTB+ build (name
# collisions and processed files end in wrong places for them to
# build):
set(BUILD_TESTING FALSE)
add_subdirectory(origin)
if(WITH_MPI)
  add_dependencies(negf mpifx)
endif()

list(APPEND PKG_CONFIG_LIBS negf)
set(PKG_CONFIG_LIBS "${PKG_CONFIG_LIBS}" PARENT_SCOPE)
