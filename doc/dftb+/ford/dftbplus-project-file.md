preprocessor:
        ../../../external/fypp/bin/fypp
macro:
        DEBUG=2
        WITH_ARPACK=1
        WITH_DFTD3=1
exclude_dir:
         ../../../prog/dftb+/lib_io/
         ../../../prog/dftb+/lib_mixer/
         ../../../prog/dftb+/lib_timedep/
         ../../../prog/dftb+/lib_derivs/
         ../../../prog/dftb+/lib_dftb/
         ../../../prog/dftb+/lib_md/
         ../../../prog/dftb+/lib_type/
exclude_file:
         ../../../prog/dftb+/lib_math/eigensolver.F90
         ../../../prog/dftb+/lib_math/factorial.F90
         ../../../prog/dftb+/lib_math/hermite.F90
         ../../../prog/dftb+/lib_math/interpolation.F90
         ../../../prog/dftb+/lib_math/lapackroutines.F90
         ../../../prog/dftb+/lib_math/make.deps
         ../../../prog/dftb+/lib_math/qm.F90
         ../../../prog/dftb+/lib_math/ranlux.F90
         ../../../prog/dftb+/lib_math/simplealgebra.F90
         ../../../prog/dftb+/prg_dftb/initprogram.F90
         ../../../prog/dftb+/prg_dftb/mainio.F90
         ../../../prog/dftb+/prg_dftb/oldcompat.F90
         ../../../prog/dftb+/prg_dftb/eigenvects.F90
         ../../../prog/dftb+/prg_dftb/inputdata.F90
         ../../../prog/dftb+/prg_dftb/parser.F90
preprocess:
        true
src_dir:
        ../../../prog/dftb+/
output_dir: ./doc
project_github: https://github.com/dftbplus/dftbplus
project_website: http://www.dftbplus.org
summary: The DFTB+ package for fast quantum mechanical atomistic simulations
author: The DFTB+ developers group
fpp_extensions:
         F90
include:
         ../../../prog/dftb+/includes/
predocmark: >
display: public
         protected
         private
proc_internals:
        false
source: true
graph: true
search: false
macro: TEST
       LOGIC=.true.
license: by
warn: true
