set(CURDIR "lib_extlibs")

set(SOURCES-FPP
  ${CURDIR}/blas.F90
  ${CURDIR}/lapack.F90
  ${CURDIR}/mpifx.F90
  ${CURDIR}/scalapackfx.F90
  ${CURDIR}/xmlf90.F90)

if(WITH_ARPACK)
  set(SOURCES-FPP ${SOURCES-FPP}
    ${CURDIR}/arpack.F90)
endif(WITH_ARPACK)

if(WITH_DFTD3)
  set(SOURCES-FPP ${SOURCES-FPP}
    ${CURDIR}/dftd3.F90)
endif(WITH_DFTD3)

if(WITH_SOCKETS)
  set(SOURCES-FPP ${SOURCES-FPP}
    ${CURDIR}/fsockets.F90)
endif(WITH_SOCKETS)

if(WITH_MPI)
  set(SOURCES-FPP ${SOURCES-FPP}
    ${CURDIR}/mpifx.f90
    ${CURDIR}/scalapackfx.f90)
endif()

if(WITH_TRANSPORT)
  set(SOURCES-FPP ${SOURCES-FPP}
    ${CURDIR}/libnegf_vars.F90
    ${CURDIR}/mat_conv.F90
    ${CURDIR}/negf_int.F90
    ${CURDIR}/poisson_int.F90)
endif(WITH_TRANSPORT)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${SOURCES-FPP} PARENT_SCOPE)
