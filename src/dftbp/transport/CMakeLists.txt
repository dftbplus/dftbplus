set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp)

if(WITH_TRANSPORT)
  list(APPEND sources-fpp
    ${curdir}/matconv.F90
    ${curdir}/negfint.F90
    ${curdir}/negfvars.F90)
endif()

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
