set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/plugin.F90)

list(APPEND sources-c
  ${curdir}/plugin.c)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
set(ALL-SOURCES-C ${ALL-SOURCES-C} ${sources-c} PARENT_SCOPE)

if(INSTALL_INCLUDE_FILES)
  install(
    FILES "${CMAKE_CURRENT_SOURCE_DIR}/dftbplus_plugin.h"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )
endif()
