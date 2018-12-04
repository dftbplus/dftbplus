set(CURDIR "lib_io")

set(SOURCES-FPP
  ${CURDIR}/charmanip.F90
  ${CURDIR}/fileid.F90
  ${CURDIR}/formatout.F90
  ${CURDIR}/hsdparser.F90
  ${CURDIR}/hsdutils.F90
  ${CURDIR}/hsdutils2.F90
  ${CURDIR}/intrinsicpr.F90
  ${CURDIR}/io.F90
  ${CURDIR}/logger.F90
  ${CURDIR}/message.F90
  ${CURDIR}/taggedoutput.F90
  ${CURDIR}/tokenreader.F90
  ${CURDIR}/xmlutils.F90)

if(WITH_SOCKETS)
  set(SOURCES-FPP ${SOURCES-FPP}
    ${CURDIR}/ipisocket.F90)
endif(WITH_SOCKETS)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${SOURCES-FPP} PARENT_SCOPE)
