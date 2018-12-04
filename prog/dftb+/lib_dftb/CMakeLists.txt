set(CURDIR "lib_dftb")

set(SOURCES-FPP
  ${CURDIR}/chargeconstr.F90
  ${CURDIR}/charges.F90
  ${CURDIR}/coulomb.F90
  ${CURDIR}/densitymatrix.F90
  ${CURDIR}/dftbplusu.F90
  ${CURDIR}/dispcommon.F90
  ${CURDIR}/dispersions.F90
  ${CURDIR}/dispiface.F90
  ${CURDIR}/dispslaterkirkw.F90
  ${CURDIR}/dispuff.F90
  ${CURDIR}/dispuffdata.F90
  ${CURDIR}/elecconstraints.F90
  ${CURDIR}/elstatpot.F90
  #${CURDIR}/emfields.F90
  ${CURDIR}/energies.F90
  ${CURDIR}/etemp.F90
  ${CURDIR}/externalcharges.F90
  ${CURDIR}/forces.F90
  ${CURDIR}/h5correction.F90
  ${CURDIR}/nonscc.F90
  ${CURDIR}/orbitalequiv.F90
  ${CURDIR}/periodic.F90
  ${CURDIR}/pmlocalisation.F90
  ${CURDIR}/populations.F90
  ${CURDIR}/potentials.F90
  ${CURDIR}/repcont.F90
  ${CURDIR}/reppoly.F90
  ${CURDIR}/repspline.F90
  ${CURDIR}/repulsive.F90
  ${CURDIR}/scc.F90
  ${CURDIR}/sccinit.F90
  ${CURDIR}/shift.F90
  ${CURDIR}/shortgamma.F90
  ${CURDIR}/sk.F90
  ${CURDIR}/slakocont.F90
  ${CURDIR}/slakoeqgrid.F90
  ${CURDIR}/sparse2dense.F90
  ${CURDIR}/spin.F90
  ${CURDIR}/spinorbit.F90
  ${CURDIR}/stress.F90
  ${CURDIR}/thirdorder.F90
  ${CURDIR}/vdwdata.F90)

if(WITH_DFTD3)
  set(SOURCES-FPP ${SOURCES-FPP}
    ${CURDIR}/dispdftd3.F90)
endif(WITH_DFTD3)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${SOURCES-FPP} PARENT_SCOPE)
