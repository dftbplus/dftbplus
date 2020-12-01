set(curdir "lib_common")

set(sources-fpp
  ${curdir}/accuracy.F90
  ${curdir}/assert.F90
  ${curdir}/atomicrad.F90
  ${curdir}/atomicmass.F90
  ${curdir}/coherence.F90
  ${curdir}/constants.F90
  ${curdir}/environment.F90
  ${curdir}/fileregistry.F90
  ${curdir}/globalenv.F90
  ${curdir}/hamiltoniantypes.F90
  ${curdir}/memman.F90
  ${curdir}/optarg.F90
  ${curdir}/schedule.F90
  ${curdir}/timer.F90
  ${curdir}/timerarray.F90
  ${curdir}/unitconversion.F90
  ${curdir}/version.F90)

if(WITH_MPI)
  list(APPEND sources-fpp
    ${curdir}/blacsenv.F90
    ${curdir}/mpienv.F90)
endif()

if(WITH_GPU)
  list(APPEND sources-fpp
    ${curdir}/gpuenv.F90)
endif()


set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
