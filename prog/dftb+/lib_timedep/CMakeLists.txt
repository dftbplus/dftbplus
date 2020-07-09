set(curdir "lib_timedep")

set(sources-fpp
  ${curdir}/dynamics.F90
  ${curdir}/linresp.F90
  ${curdir}/linrespcommon.F90
  ${curdir}/linrespgrad.F90
  ${curdir}/linresptypes.F90
  ${curdir}/pprpa.F90
  ${curdir}/transcharges.F90)

if(WITH_ARPACK)
  list(APPEND sources-fpp
    ${curdir}/rs_linresp.F90)
endif(WITH_ARPACK)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
