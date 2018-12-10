set(curdir "lib_extlibs")

set(sources-fpp
  ${curdir}/blas.F90
  ${curdir}/lapack.F90
  ${curdir}/mpifx.F90
  ${curdir}/scalapackfx.F90
  ${curdir}/xmlf90.F90)

if(WITH_ARPACK)
  list(APPEND sources-fpp ${curdir}/arpack.F90)
endif()

if(WITH_DFTD3)
  list(APPEND sources-fpp ${curdir}/dftd3.F90)
endif()

if(WITH_SOCKETS)
  list(APPEND sources-fpp ${curdir}/fsockets.F90)
endif(WITH_SOCKETS)

if(WITH_MPI)
  list(APPEND sources-fpp
    ${curdir}/mpifx.f90
    ${curdir}/scalapackfx.f90)
endif()

if(WITH_TRANSPORT)
  list(APPEND sources-fpp
    ${curdir}/libnegf_vars.F90
    ${curdir}/mat_conv.F90
    ${curdir}/negf_int.F90
    ${curdir}/poisson_int.F90)
endif()

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
