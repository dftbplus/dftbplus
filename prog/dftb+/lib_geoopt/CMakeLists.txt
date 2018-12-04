set(CURDIR "lib_geoopt")

set(SOURCES-FPP
  ${CURDIR}/conjgrad.F90
  ${CURDIR}/gdiis.F90
  ${CURDIR}/geoopt.F90
  ${CURDIR}/lbfgs.F90
  ${CURDIR}/linemin.F90
  ${CURDIR}/steepdesc.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${SOURCES-FPP} PARENT_SCOPE)
