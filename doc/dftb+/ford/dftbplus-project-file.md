preprocessor:
        ../../../external/fypp/bin/fypp
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
source: false
graph: false
search: false
macro: TEST
       LOGIC=.true.
license: by
warn: true
