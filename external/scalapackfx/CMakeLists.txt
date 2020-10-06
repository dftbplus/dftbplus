set(BUILD_EXPORTED_TARGETS_ONLY TRUE)
add_subdirectory(origin)
add_library(ScalapackFx::ScalapackFx INTERFACE IMPORTED GLOBAL)
target_link_libraries(ScalapackFx::ScalapackFx INTERFACE ScalapackFx)
