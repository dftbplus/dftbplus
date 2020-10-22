macro(before_source_config)
  set(BUILD_EXPORTED_TARGETS_ONLY TRUE)
endmacro()

macro(after_source_config)
  add_library(MpiFx::MpiFx INTERFACE IMPORTED GLOBAL)
  target_link_libraries(MpiFx::MpiFx INTERFACE MpiFx)
endmacro()

set(MPIFX_GIT_REPOSITORY "https://github.com/dftbplus/mpifx.git")
set(MPIFX_GIT_TAG "1.0")

dftbp_config_hybrid_source_dependency(${MPIFX_GIT_REPOSITORY} ${MPIFX_GIT_TAG})
