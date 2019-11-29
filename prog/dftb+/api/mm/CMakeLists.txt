set(curdir "api/mm")

set(sources-fpp
  ${curdir}/capi.F90
  ${curdir}/dftbplus.F90
  ${curdir}/hsdapi.F90
  ${curdir}/mmapi.F90
  ${curdir}/qdepextpotgenc.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
