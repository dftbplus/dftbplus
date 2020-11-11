add_subdirectory(dftb+)
if(NOT BUILD_EXPORTED_TARGETS_ONLY)
  add_subdirectory(misc)
  add_subdirectory(modes)
  add_subdirectory(waveplot)

  if(WITH_TRANSPORT)
    add_subdirectory(transporttools)
    add_subdirectory(phonons)
  endif()
endif()

# Propagate PKG_CONFIG_LIBS upwards
set(PKG_CONFIG_LIBS "${PKG_CONFIG_LIBS}" PARENT_SCOPE)
