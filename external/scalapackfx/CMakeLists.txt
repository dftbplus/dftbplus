macro(before_source_config)
  set(BUILD_EXPORTED_TARGETS_ONLY TRUE)
endmacro()

macro(after_source_config)
  add_library(ScalapackFx::ScalapackFx INTERFACE IMPORTED GLOBAL)
  target_link_libraries(ScalapackFx::ScalapackFx INTERFACE ScalapackFx)
endmacro()

dftbp_config_hybrid_source_dependency(ScalapackFx SCALAPACKFX)
