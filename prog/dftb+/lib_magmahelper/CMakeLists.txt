set(curdir "lib_magmahelper")

set(sources-c
  ${curdir}/count_gpus.c)

set(sources-fpp
  ${curdir}/magma_device_info.f90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
set(ALL-SOURCES-C ${ALL-SOURCES-C} ${sources-c} PARENT_SCOPE)
