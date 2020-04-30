set(INSTALL_EXPORT_NAME dftbplus-targets)
set(MPIFX_INCLUDE_DIR "${CMAKE_BINARY_DIR}/external/mpifx/origin/lib/include" CACHE STRING
  "Where to find the MPIFX include directory")
add_subdirectory(origin)
set(EXPORTED_COMPILED_LIBRARIES ${EXPORTED_COMPILED_LIBRARIES} negf PARENT_SCOPE)
if(WITH_MPI)
  add_dependencies(negf mpifx)
endif()
