set(INSTALL_EXPORT_NAME dftbplus-targets)
set(MPIFX_INCLUDE_DIR "${CMAKE_BINARY_DIR}/external/mpifx/origin/lib/include" CACHE STRING
  "Where to find the MPIFX include directory")
# prevent libNEGF testers being built if inside DFTB+ build (name
# collisions and processed files end in wrong places for them to
# build):
set(BUILD_TESTING CACHE "Whether the tests should be built" FALSE)
add_subdirectory(origin)
set(EXPORTED_COMPILED_LIBRARIES ${EXPORTED_COMPILED_LIBRARIES} negf PARENT_SCOPE)
if(WITH_MPI)
  add_dependencies(negf mpifx)
endif()
