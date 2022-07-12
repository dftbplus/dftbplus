set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/dynamicsrestart.F90
  ${curdir}/linresp.F90
  ${curdir}/linrespcommon.F90
  ${curdir}/linrespgrad.F90
  ${curdir}/linresptypes.F90
  ${curdir}/pprpa.F90
  ${curdir}/transcharges.F90
  ${curdir}/timeprop.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
