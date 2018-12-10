set(curdir "prg_dftb")

set(sources-fpp
  ${curdir}/dftbplus.F90
  ${curdir}/eigenvects.F90
  ${curdir}/initprogram.F90
  ${curdir}/inputconversion.F90
  ${curdir}/inputdata.F90
  ${curdir}/main.F90
  ${curdir}/mainio.F90
  ${curdir}/oldcompat.F90
  ${curdir}/parser.F90
  ${curdir}/solvertypes.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
