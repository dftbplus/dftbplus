set(curdir "lib_timedep")

set(sources-fpp
  ${curdir}/linresp.F90
  ${curdir}/linrespcommon.F90
  ${curdir}/linrespgrad.F90
  ${curdir}/transcharges.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
