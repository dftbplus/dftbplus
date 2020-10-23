macro(before_source_config)
  set(ENABLE_SCALAPACK_MPI ${WITH_MPI})
  set(ENABLE_ELSI ${WITH_ELSI})
  set(ENABLE_C_API OFF)
  set(SCALAPACK_LIBRARIES "${SCALAPACK_LIBRARY}")
  set(MKL_LIBRARIES "")
  set(BUILD_TESTING FALSE)
endmacro()

macro(after_source_config)
  add_library(Mbd::Mbd INTERFACE IMPORTED GLOBAL)
  target_link_libraries(Mbd::Mbd INTERFACE Mbd)
endmacro()

set(MBD_GIT_REPOSITORY "https://github.com/jhrmnn/libmbd.git")
# Use utils/test/check_submodule_commits to update commit hash below
set(MBD_GIT_TAG "3e70b0cc64557b1bf06467bb2d2d50cee8e2e7ef")

dftbp_config_hybrid_source_dependency(${MBD_GIT_REPOSITORY} ${MBD_GIT_TAG})
