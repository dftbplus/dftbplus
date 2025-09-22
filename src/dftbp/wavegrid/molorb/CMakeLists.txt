set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/molorb.F90
  ${curdir}/parallel.F90
  ${curdir}/offloaded.F90
  ${curdir}/spharmonics.F90
  ${curdir}/types.F90
)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
