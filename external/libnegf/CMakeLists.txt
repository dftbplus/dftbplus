macro(before_source_config)
  set(BUILD_TESTING FALSE)
endmacro()

macro(after_source_config)
  add_library(Negf::Negf INTERFACE IMPORTED GLOBAL)
  target_link_libraries(Negf::Negf INTERFACE Negf)
endmacro()

set(NEGF_GIT_REPOSITORY "https://github.com/libnegf/libnegf.git")
# Use utils/test/check_submodule_commits to update commit hash below
set(NEGF_GIT_TAG "3b9a727868182002d7570d2676aaf9bccb299bcc")

dftbp_config_hybrid_source_dependency(${NEGF_GIT_REPOSITORY} ${NEGF_GIT_TAG})
