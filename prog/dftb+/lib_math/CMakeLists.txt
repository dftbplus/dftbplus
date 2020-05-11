set(curdir "lib_math")

set(sources-fpp
  ${curdir}/angmomentum.F90
  ${curdir}/bisect.F90
  ${curdir}/blasroutines.F90
  ${curdir}/degeneracy.F90
  ${curdir}/determinant.F90
  ${curdir}/eigensolver.F90
  ${curdir}/erfcalc.F90
  ${curdir}/errorfunction.F90
  ${curdir}/f08math.F90
  ${curdir}/factorial.F90
  ${curdir}/hermite.F90
  ${curdir}/interpolation.F90
  ${curdir}/lapackroutines.F90
  ${curdir}/qm.F90
  ${curdir}/quaternions.F90
  ${curdir}/randomgenpool.F90
  ${curdir}/ranlux.F90
  ${curdir}/simplealgebra.F90
  ${curdir}/sorting.F90)

if(WITH_MPI)
  list(APPEND sources-fpp
    ${curdir}/scalafxext.F90)
endif()

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
