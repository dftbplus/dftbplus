set(curdir "lib_io")

set(sources-fpp
  ${curdir}/charmanip.F90
  ${curdir}/fileid.F90
  ${curdir}/formatout.F90
  ${curdir}/hsdparser.F90
  ${curdir}/hsdutils.F90
  ${curdir}/hsdutils2.F90
  ${curdir}/intrinsicpr.F90
  ${curdir}/logger.F90
  ${curdir}/message.F90
  ${curdir}/taggedoutput.F90
  ${curdir}/tokenreader.F90
  ${curdir}/xmlutils.F90)

if(WITH_SOCKETS)
  list(APPEND sources-fpp ${curdir}/ipisocket.F90)
endif()

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
