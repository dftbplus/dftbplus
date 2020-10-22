macro(before_source_config)
  set(BUILD_EXPORTED_TARGETS_ONLY TRUE)
endmacro()

macro(after_source_config)
  add_library(MpiFx::MpiFx INTERFACE IMPORTED GLOBAL)
  target_link_libraries(MpiFx::MpiFx INTERFACE MpiFx)
endmacro()

set(MPIFX_GIT_REPOSITORY "https://github.com/dftbplus/mpifx.git")
# Use utils/test/check_submodule_commits to update commit hash below
set(MPIFX_GIT_TAG "0cb07ee08cbb20f3f7bb2527152a4ec317c579ad")

dftbp_config_hybrid_source_dependency(${MPIFX_GIT_REPOSITORY} ${MPIFX_GIT_TAG})
