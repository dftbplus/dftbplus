preprocessor:
        ../../../external/fypp/bin/fypp
macro:
        DEBUG=2
	WITH_MPI=0
        WITH_ARPACK=0
        WITH_DFTD3=1
        WITH_SOCKETS=0
preprocess:
        true
src_dir:
        ../../../api/mm/
	../../../test/api/mm
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
proc_internals:
        false
source: true
graph: true
search: false
macro: TEST
       LOGIC=.true.
license: by
warn: true
