set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/andersonmixer.F90
  ${curdir}/broydenmixer.F90
  ${curdir}/diismixer.F90
  ${curdir}/mixer.F90
  ${curdir}/simplemixer.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
