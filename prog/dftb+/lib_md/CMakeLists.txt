set(CURDIR "lib_md")

set(SOURCES-FPP
  ${CURDIR}/andersentherm.F90
  ${CURDIR}/berendsentherm.F90
  ${CURDIR}/dummytherm.F90
  ${CURDIR}/extlagrangian.F90
  ${CURDIR}/mdcommon.F90
  ${CURDIR}/mdintegrator.F90
  ${CURDIR}/nhctherm.F90
  ${CURDIR}/tempprofile.F90
  ${CURDIR}/thermostat.F90
  ${CURDIR}/velocityverlet.F90
  ${CURDIR}/xlbomd.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${SOURCES-FPP} PARENT_SCOPE)
