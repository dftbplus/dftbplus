set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/arpack.F90
  ${curdir}/blas.F90
  ${curdir}/chimes.F90
  ${curdir}/ddcosmo.F90
  ${curdir}/dftd4refs.F90
  ${curdir}/elpa.F90
  ${curdir}/elsiiface.F90
  ${curdir}/lapack.F90
  ${curdir}/lebedev.F90
  ${curdir}/magma.F90
  ${curdir}/mpifx.F90
  ${curdir}/plumed.F90
  ${curdir}/poisson.F90
  ${curdir}/scalapackfx.F90
  ${curdir}/sdftd3.F90
  ${curdir}/tblite.F90
  ${curdir}/xmlf90.F90)

if(WITH_SOCKETS)
  list(APPEND sources-fpp ${curdir}/fsockets.F90)
endif(WITH_SOCKETS)

if(WITH_TRANSPORT)
  list(APPEND sources-fpp
    ${curdir}/negf.F90)
endif()

if(WITH_MAGMA)
  list(APPEND sources-c
    ${curdir}/magmac.c)
endif()

if(WITH_MBD)
  list(APPEND sources-fpp ${curdir}/mbd.F90)
endif()

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
set(ALL-SOURCES-C ${ALL-SOURCES-C} ${sources-c} PARENT_SCOPE)
