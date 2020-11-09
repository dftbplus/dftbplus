set(curdir "lib_geoopt")

set(sources-fpp
  ${curdir}/conjgrad.F90
  ${curdir}/fire.F90
  ${curdir}/gdiis.F90
  ${curdir}/geoopt.F90
  ${curdir}/lbfgs.F90
  ${curdir}/linemin.F90
  ${curdir}/steepdesc.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
