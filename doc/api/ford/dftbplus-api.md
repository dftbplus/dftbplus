macro:
        WITH_MPI=0
        WITH_ARPACK=0
        WITH_DFTD3=1
        WITH_SOCKETS=0
        RELEASE=20.2
preprocess: true
src_dir:
        ../../../prog/dftb+/api/mm
        ../../../test/api/mm
output_dir: ./doc
project_github: https://github.com/dftbplus/dftbplus
project_website: http://www.dftbplus.org
summary: The DFTB+ package for fast quantum mechanical atomistic simulations
author: The DFTB+ developers group
preprocessor: ../../../external/fypp/bin/fypp
include: ../../../prog/dftb+/include
predocmark: >
display: public
         protected
proc_internals:
        false
source: true
graph: true
search: false
license: by
warn: true
