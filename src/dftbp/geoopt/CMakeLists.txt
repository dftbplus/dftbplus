set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/conjgrad.F90
  ${curdir}/filter.F90
  ${curdir}/fire.F90
  ${curdir}/gdiis.F90
  ${curdir}/geoopt.F90
  ${curdir}/lbfgs.F90
  ${curdir}/lbfgs2.F90
  ${curdir}/linemin.F90
  ${curdir}/optimizer.F90
  ${curdir}/package.F90
  ${curdir}/rationalfunc.F90
  ${curdir}/steepdesc.F90)

add_subdirectory(deprecated)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
