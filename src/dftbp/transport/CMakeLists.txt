set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp)

if(WITH_TRANSPORT)
  list(APPEND sources-fpp
    ${curdir}/matconv.F90
    ${curdir}/negfint.F90)
endif()
list(APPEND sources-fpp
  ${curdir}/negfvars.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
