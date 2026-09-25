set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/basis.F90
  ${curdir}/slater.F90
  ${curdir}/lut.F90
  ${curdir}/spharmonics.F90
  ${curdir}/orbital.F90
)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
