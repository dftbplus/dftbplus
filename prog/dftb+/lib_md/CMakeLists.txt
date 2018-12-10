set(curdir "lib_md")

set(sources-fpp
  ${curdir}/andersentherm.F90
  ${curdir}/berendsentherm.F90
  ${curdir}/dummytherm.F90
  ${curdir}/extlagrangian.F90
  ${curdir}/mdcommon.F90
  ${curdir}/mdintegrator.F90
  ${curdir}/nhctherm.F90
  ${curdir}/tempprofile.F90
  ${curdir}/thermostat.F90
  ${curdir}/velocityverlet.F90
  ${curdir}/xlbomd.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
