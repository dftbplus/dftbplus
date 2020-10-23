macro(before_source_config)
  set(BUILD_EXPORTED_TARGETS_ONLY TRUE)
endmacro()

macro(after_source_config)
  add_library(ScalapackFx::ScalapackFx INTERFACE IMPORTED GLOBAL)
  target_link_libraries(ScalapackFx::ScalapackFx INTERFACE ScalapackFx)
endmacro()

set(SCALAPACKFX_GIT_REPOSITORY "https://github.com/dftbplus/scalapackfx.git")
# Use utils/test/check_submodule_commits to update commit hash below
set(SCALAPACKFX_GIT_TAG "cd758a283a89905c5a4fd3e4a42305c287559bae")

dftbp_config_hybrid_source_dependency(${SCALAPACKFX_GIT_REPOSITORY} ${SCALAPACKFX_GIT_TAG})
