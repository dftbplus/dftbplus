set(CURDIR "prg_dftb")

set(SOURCES-FPP
  ${CURDIR}/dftbplus.F90
  ${CURDIR}/eigenvects.F90
  ${CURDIR}/initprogram.F90
  ${CURDIR}/inputconversion.F90
  ${CURDIR}/inputdata.F90
  ${CURDIR}/main.F90
  ${CURDIR}/mainio.F90
  ${CURDIR}/oldcompat.F90
  ${CURDIR}/parser.F90
  ${CURDIR}/solvertypes.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${SOURCES-FPP} PARENT_SCOPE)
