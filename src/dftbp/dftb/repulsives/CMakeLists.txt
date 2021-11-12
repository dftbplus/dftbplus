set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

list(APPEND ALL-SOURCES-FPP
  ${curdir}/pairrepulsive.F90
  ${curdir}/polyrep.F90
  ${curdir}/repulsive.F90
  ${curdir}/splinerep.F90
  ${curdir}/twobodyrep.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} PARENT_SCOPE)
