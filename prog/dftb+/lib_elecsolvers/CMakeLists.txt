set(curdir "lib_elecsolvers")

set(sources-fpp
  ${curdir}/elecsolvers.F90
  ${curdir}/elecsolvertypes.F90
  ${curdir}/elsicsc.F90
  ${curdir}/elsisolver.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
