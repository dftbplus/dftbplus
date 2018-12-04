set(CURDIR "lib_math")

set(SOURCES-FPP
  ${CURDIR}/angmomentum.F90
  ${CURDIR}/bisect.F90
  ${CURDIR}/blasroutines.F90
  ${CURDIR}/determinant.F90
  ${CURDIR}/eigensolver.F90
  ${CURDIR}/erfcalc.F90
  ${CURDIR}/errorfunction.F90
  ${CURDIR}/factorial.F90
  ${CURDIR}/hermite.F90
  ${CURDIR}/interpolation.F90
  ${CURDIR}/lapackroutines.F90
  ${CURDIR}/qm.F90
  ${CURDIR}/randomgenpool.F90
  ${CURDIR}/ranlux.F90
  ${CURDIR}/simplealgebra.F90
  ${CURDIR}/sorting.F90)

if(WITH_MPI)
  set(SOURCES-FPP ${SOURCES-FPP}
    ${CURDIR}/scalafxext.F90)
endif(WITH_MPI)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${SOURCES-FPP} PARENT_SCOPE)
