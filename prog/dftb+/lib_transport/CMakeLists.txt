set(curdir "lib_transport")

set(sources-fpp)

if(WITH_TRANSPORT)
  list(APPEND sources-fpp
    ${curdir}/negfint.F90
    ${curdir}/negfvars.F90
    ${curdir}/matconv.F90)
endif()

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
