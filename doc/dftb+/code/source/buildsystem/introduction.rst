DFTB+ build system
==================

Cmake
-----

The main (Fortran based) programs of the DFTB+ project use cmake to
organize compilation of the project.

There are project specific customisations inside the `cmake/` folder::


  ├── cmake/
  │   ├── DftbPlusUtils.cmake
  │   └── Modules/
  │       ├── CustomLibraryFinder.cmake
  │       ├
  │
  ├── external/
  │   ├── xmlf90/
  │   ├


  * `DftbPlusUtils.cmake` contains utility functions and also
    definitions for `-D` cmake flags.

  * `Modules` currently contains finder routines for various optional
    and required library packages (most of the option libraries are
    currently located in the `external` sub-directory).

Pre-processing of source code
-----------------------------

Fortran source code with the `.F90` extension is pre-processed into
`.f90` files before compilation. This is done using the `fypp
<https://github.com/aradi/fypp>`_ Python-based preprocessor (see also
it's `documentation <https://fypp.readthedocs.io/>`_

To retain intermediate `.f90` files, select the appropriate option for
cmake ::

  cmake --debug-trycompile -B_build .

Then, after running `make` the preprocesssed intermediate files will
be located inside `_build/src/dftbp/`

Use of fypp
-----------

See the detailed `fypp documentation
<https://fypp.readthedocs.io/en/stable/>`_ for details, but in brief,
fypp provices Python-based metaprogramming capabilities for Fortran,
allowing code macros, variable substitutions, looped and conditional
generation of code, etc.

Variable conventions
~~~~~~~~~~~~~~~~~~~~

Fortran-like choices should be preferred, for example if defining
arguments for a generated routine, choose naming conventions which are
familar to Fortran programmers ::

  #:for TYPE, KIND VAR in [("real", "dp", "A"), ("real", "sp", "B"), ("complex", "dp", "C")]
  
  ${TYPE}$(${KIND}$) :: variable${VAR}$

  #:endfor

would generate code that looks like ::

  real(dp) :: variableA
  real(sp) :: variableB
  complex(dp) :: variableC


Suggestions for project submodules
----------------------------------

External projects used by DFTB+ are located in the `external/`
folder. New projects should usually be added as git `submodules`::

  external/
  ├── projectname
  │   └── origin

You should add a CMakeLists.txt file in the `projectname/` directory
of the module to control its compilation.

Then customize the top level `CMakeLists.txt` to find the library, add
a link dependency in `src/dftbp/CMakeLists.txt` and appropriate
wrappers in `src/dftbp/extlibs/` and a dependency requirement to
`utils/export/dftbplus-config.cmake.in`.

The `utils/test/check_submodule_commits` script can be used to
check/update the git tag in the `CMakeLists.txt` and for consistency
between the provided url and the local `.git/config` address.

Note that the url used in `CMakeLists.txt` should point to an https
address for the submodule's repository, not a git address for the
regression testing pipeline to be able to download them.
