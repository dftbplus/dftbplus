set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/dmsolvertypes.F90
  ${curdir}/elecsolvers.F90
  ${curdir}/elecsolvertypes.F90
  ${curdir}/elsicsc.F90
  ${curdir}/elsisolver.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
