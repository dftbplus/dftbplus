# This CMakeLists.txt should be invoked via dftbp_config_hybrid_dependency()!

set(BUILD_EXPORTED_TARGETS_ONLY TRUE)

add_subdirectory(${SCALAPACKFX_SOURCE_DIR} ${SCALAPACKFX_BINARY_DIR})

add_library(ScalapackFx::ScalapackFx INTERFACE IMPORTED GLOBAL)
target_link_libraries(ScalapackFx::ScalapackFx INTERFACE ScalapackFx)
