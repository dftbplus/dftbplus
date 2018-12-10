set(curdir "lib_type")

set(sources-fpp
  ${curdir}/commontypes.F90
  ${curdir}/densedescr.F90
  ${curdir}/dynneighlist.F90
  ${curdir}/latpointiter.F90
  ${curdir}/linkedlist.F90
  ${curdir}/linkedlisti0.F90
  ${curdir}/linkedlisti1.F90
  ${curdir}/linkedlistlc0.F90
  ${curdir}/linkedlistmc0.F90
  ${curdir}/linkedlistr0.F90
  ${curdir}/linkedlistr1.F90
  ${curdir}/linkedlistr2.F90
  ${curdir}/linkedlists0.F90
  ${curdir}/oldskdata.F90
  ${curdir}/orbitals.F90
  ${curdir}/parallelks.F90
  ${curdir}/typegeometry.F90
  ${curdir}/typegeometryhsd.F90
  ${curdir}/wrapped_intrinsics.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
