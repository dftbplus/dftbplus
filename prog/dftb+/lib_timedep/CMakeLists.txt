set(CURDIR "lib_timedep")

set(SOURCES-FPP
  ${CURDIR}/linresp.F90)

if(WITH_ARPACK)
  set(SOURCES-FPP ${SOURCES-FPP}
    ${CURDIR}/linrespcommon.F90
    ${CURDIR}/linrespgrad.F90)
endif(WITH_ARPACK)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${SOURCES-FPP} PARENT_SCOPE)
