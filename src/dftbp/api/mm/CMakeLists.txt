set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/capi.F90
  ${curdir}/dftbplus.F90
  ${curdir}/hsdapi.F90
  ${curdir}/mmapi.F90
  ${curdir}/qdepextpotgenc.F90)

if(NOT WITH_EXTERNALHAM)
  list(APPEND sources-fpp
    ${curdir}/dummyextham.F90
    )
endif()

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
