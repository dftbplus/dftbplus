set(curdir "lib_timedep")

set(sources-fpp
  ${curdir}/linresp.F90
  ${curdir}/pprpa.F90)

if(WITH_ARPACK)
  list(APPEND sources-fpp ${curdir}/linrespgrad.F90)
endif(WITH_ARPACK)

list(APPEND sources-fpp
  ${curdir}/linrespcommon.F90
  ${curdir}/transcharges.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
