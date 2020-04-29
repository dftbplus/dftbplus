set(curdir "lib_solvation")

set(sources-fpp
  ${curdir}/born.F90
  ${curdir}/borndata.F90
  ${curdir}/solventdata.F90
  ${curdir}/solvinput.F90
  ${curdir}/solvation.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)

