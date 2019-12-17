set(curdir "lib_reks")

set(sources-fpp
  ${curdir}/rekscommon.F90
  ${curdir}/reksen.F90
  ${curdir}/reksfon.F90
  ${curdir}/reksproperty.F90
  ${curdir}/reksvar.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
