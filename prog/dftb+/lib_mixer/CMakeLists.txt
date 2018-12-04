set(CURDIR "lib_mixer")

set(SOURCES-FPP
  ${CURDIR}/andersonmixer.F90
  ${CURDIR}/broydenmixer.F90
  ${CURDIR}/diismixer.F90
  ${CURDIR}/mixer.F90
  ${CURDIR}/simplemixer.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${SOURCES-FPP} PARENT_SCOPE)
