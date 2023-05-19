set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/accuracy.F90
  ${curdir}/assert.F90
  ${curdir}/atomicrad.F90
  ${curdir}/atomicmass.F90
  ${curdir}/coherence.F90
  ${curdir}/constants.F90
  ${curdir}/envcheck.F90
  ${curdir}/environment.F90
  ${curdir}/file.F90
  ${curdir}/filesystem.F90
  ${curdir}/globalenv.F90
  ${curdir}/hamiltoniantypes.F90
  ${curdir}/memman.F90
  ${curdir}/optarg.F90
  ${curdir}/release.F90
  ${curdir}/parallel.F90
  ${curdir}/schedule.F90
  ${curdir}/status.F90
  ${curdir}/timer.F90
  ${curdir}/timerarray.F90
  ${curdir}/unitconversion.F90
  ${curdir}/version.F90)

list(APPEND sources-c
  ${curdir}/envcheck.c)

if(WITH_MPI)
  list(APPEND sources-fpp
    ${curdir}/blacsenv.F90
    ${curdir}/mpienv.F90)
endif()

if(WITH_MAGMA)
  list(APPEND sources-fpp
    ${curdir}/gpuenv.F90)
endif()


set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
set(ALL-SOURCES-C ${ALL-SOURCES-C} ${sources-c} PARENT_SCOPE)
