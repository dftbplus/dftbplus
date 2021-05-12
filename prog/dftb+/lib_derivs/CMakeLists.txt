set(curdir "lib_derivs")

set(sources-fpp
  ${curdir}/fermihelper.F90
  ${curdir}/numderivs2.F90
  ${curdir}/rotatedegen.F90
  ${curdir}/staticperturb.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
