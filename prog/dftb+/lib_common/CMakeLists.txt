set(CURDIR "lib_common")

set(SOURCES-FPP
  ${CURDIR}/accuracy.F90
  ${CURDIR}/assert.F90
  ${CURDIR}/constants.F90
  ${CURDIR}/environment.F90
  ${CURDIR}/fileregistry.F90
  ${CURDIR}/globalenv.F90
  ${CURDIR}/memman.F90
  ${CURDIR}/optarg.F90
  ${CURDIR}/schedule.F90
  ${CURDIR}/timer.F90
  ${CURDIR}/timerarray.F90
  ${CURDIR}/unitconversion.F90)

if(WITH_MPI)
  set(SOURCES-FPP ${SOURCES-FPP}
    ${CURDIR}/blacsenv.F90
    ${CURDIR}/mpienv.F90)
endif(WITH_MPI)


set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${SOURCES-FPP} PARENT_SCOPE)
