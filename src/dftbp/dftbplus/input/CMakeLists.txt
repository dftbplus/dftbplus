set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

list(APPEND ALL-SOURCES-FPP
  ${curdir}/fileaccess.F90
  ${curdir}/geoopt.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} PARENT_SCOPE)
