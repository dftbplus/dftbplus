src_dir: ../../../prog/dftb+/
exclude_dir:
        ../../../prog/dftb+/lib_dftb/
        ../../../prog/dftb+/lib_io/
        ../../../prog/dftb+/lib_mixer/
        ../../../prog/dftb+/lib_extlibs/
        ../../../prog/dftb+/lib_timedep/
        ../../../prog/dftb+/lib_derivs/
        ../../../prog/dftb+/lib_md/
        ../../../prog/dftb+/lib_type/
        ../../../prog/dftb+/prg_dftb/
        ../../../prog/dftb+/lib_math/
exclude_file:
../../../prog/dftb+/lib_common/constants.F90
../../../prog/dftb+/lib_common/memman.F90
../../../prog/dftb+/lib_common/unitconversion.F90
../../../prog/dftb+/lib_common/assert.F90
../../../prog/dftb+/lib_common/optarg.F90
output_dir: ./doc
project_github: https://github.com/dftbplus/dftbplus
project_website: http://www.dftbplus.org
summary: The DFTB+ package for fast quantum mechanical atomistic simulations
author: The DFTB+ developers group
fpp_extensions:
         F90
preprocessor: ../../../external/fypp/bin/fypp
preprocess: true

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
search: true
macro: TEST
       LOGIC=.true.
license: by
warn: true
