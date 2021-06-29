set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/born.F90
  ${curdir}/cm5.F90
  ${curdir}/cosmo.F90
  ${curdir}/gbsafile.F90
  ${curdir}/sasa.F90
  ${curdir}/solvation.F90
  ${curdir}/solvdata.F90
  ${curdir}/solventdata.F90
  ${curdir}/solvinput.F90
  ${curdir}/solvparser.F90
)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)

