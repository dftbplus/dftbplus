set(CURDIR "lib_type")

set(SOURCES-FPP
  ${CURDIR}/commontypes.F90
  ${CURDIR}/densedescr.F90
  ${CURDIR}/dynneighlist.F90
  ${CURDIR}/latpointiter.F90
  ${CURDIR}/linkedlist.F90
  ${CURDIR}/linkedlisti0.F90
  ${CURDIR}/linkedlisti1.F90
  ${CURDIR}/linkedlistlc0.F90
  ${CURDIR}/linkedlistmc0.F90
  ${CURDIR}/linkedlistr0.F90
  ${CURDIR}/linkedlistr1.F90
  ${CURDIR}/linkedlistr2.F90
  ${CURDIR}/linkedlists0.F90
  ${CURDIR}/oldskdata.F90
  ${CURDIR}/orbitals.F90
  ${CURDIR}/parallelks.F90
  ${CURDIR}/typegeometry.F90
  ${CURDIR}/typegeometryhsd.F90
  ${CURDIR}/wrapped_intrinsics.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${SOURCES-FPP} PARENT_SCOPE)
