macro:
        DEBUG=0
        WITH_MPI=1
        WITH_ARPACK=1
        WITH_DFTD3=1
        WITH_SOCKETS=1
        RELEASE=22.1
preprocess: true
src_dir:
        ../../../app/dftb+
        ../../../app/modes
        ../../../app/waveplot
        ../../../app/misc
output_dir: ./doc
project_github: https://github.com/dftbplus/dftbplus
project_website: http://www.dftbplus.org
summary: The DFTB+ package for fast quantum mechanical atomistic simulations
author: The DFTB+ developers group
preprocessor: ../../../external/fypp/bin/fypp
include: ../../../src/dftbp/include
         ../../../external/dftd4refs
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
