set(curdir "lib_derivs")

set(sources-fpp
  ${curdir}/finiteT.F90
  ${curdir}/numderivs2.F90
  ${curdir}/perturbationDerivs.F90
  ${curdir}/perturbationk.F90
  ${curdir}/perturbationx.F90
  ${curdir}/rotateDegenerateOrbitals.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
