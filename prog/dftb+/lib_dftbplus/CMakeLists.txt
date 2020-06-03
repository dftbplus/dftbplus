set(curdir "lib_dftbplus")

set(sources-fpp
  ${curdir}/eigenvects.F90
  ${curdir}/elstattypes.F90
  ${curdir}/forcetypes.F90
  ${curdir}/gpuinfo.F90
  ${curdir}/hsdhelpers.F90
  ${curdir}/initprogram.F90
  ${curdir}/inputconversion.F90
  ${curdir}/inputdata.F90
  ${curdir}/main.F90
  ${curdir}/mainapi.F90
  ${curdir}/mainio.F90
  ${curdir}/oldcompat.F90
  ${curdir}/parser.F90
  ${curdir}/qdepextpotgen.F90
  ${curdir}/qdepextpotproxy.F90
  ${curdir}/specieslist.F90
  ${curdir}/transportio.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
