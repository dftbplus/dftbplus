macro(before_source_config)
  set(BUILD_EXPORTED_TARGETS_ONLY TRUE)
endmacro()

macro(after_source_config)
  add_library(MpiFx::MpiFx INTERFACE IMPORTED GLOBAL)
  target_link_libraries(MpiFx::MpiFx INTERFACE MpiFx)
endmacro()

dftbp_config_hybrid_source_dependency(MpiFx MPIFX)
