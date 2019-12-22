set(sources-f77
  d1mach.f
  d9gmic.f
  d9gmit.f
  d9lgic.f
  d9lgit.f
  d9lgmc.f
  dcsevl.f
  dgamic.f
  dgamit.f
  dgamlm.f
  dgamln.f
  dgamma.f
  dgamr.f
  dlgams.f
  dlngam.f
  fdump.f
  i1mach.f
  initds.f
  j4save.f
  xerclr.f
  xercnt.f
  xerhlt.f
  xermsg.f
  xerprn.f
  xersve.f
  xgetf.f
  xgetua.f
  xsetf.f
  zabs.f
  zacai.f
  zairy.f
  zasyi.f
  zbesi.f
  zbinu.f
  zbknu.f
  zbuni.f
  zdiv.f
  zexp.f
  zkscl.f
  zlog.f
  zmlri.f
  zmlt.f
  zrati.f
  zs1s2.f
  zseri.f
  zshch.f
  zsqrt.f
  zuchk.f
  zunhj.f
  zuni1.f
  zuni2.f
  zunik.f
  zuoik.f
  zwrsk.f)

add_library(slatec ${sources-f77})

install(TARGETS slatec DESTINATION ${INSTALL_LIB_DIR} EXPORT dftbplus-targets)

set(EXPORTED_COMPILED_LIBRARIES ${EXPORTED_COMPILED_LIBRARIES} slatec PARENT_SCOPE)
