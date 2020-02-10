set(curdir "lib_xtb")

set(sources-fpp
  ${curdir}/coordinationnumber.F90
  ${curdir}/gfn1param.F90
  ${curdir}/gtocont.F90
  ${curdir}/gtoints.F90
  ${curdir}/repsimple.F90
  ${curdir}/slater.F90
  ${curdir}/xtbh0.F90
  ${curdir}/xtbinput.F90
  ${curdir}/xtbparam.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
