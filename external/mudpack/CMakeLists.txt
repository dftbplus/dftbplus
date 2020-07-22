set(sources-f77
  dmud3.f
  dmuh3.f
  dmud3ln.f
  dmud3pn.f
  dmud3sp.f
  dmudcom.f)

add_library(mudpack ${sources-f77})

if(WITH_API OR BUILD_SHARED_LIBS)
  install(TARGETS mudpack DESTINATION ${INSTALL_LIB_DIR} EXPORT dftbplus-targets)
endif()

set(EXPORTED_COMPILED_LIBRARIES ${EXPORTED_COMPILED_LIBRARIES} mudpack PARENT_SCOPE)
