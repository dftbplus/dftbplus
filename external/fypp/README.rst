*********************************************
Fypp — Python powered Fortran metaprogramming
*********************************************

.. image:: https://travis-ci.org/aradi/fypp.svg?branch=develop
           :target: https://travis-ci.org/aradi/fypp

Fypp is a Python powered preprocessor. It can be used for any programming
languages but its primary aim is to offer a Fortran preprocessor, which helps to
extend Fortran with condititional compiling and template metaprogramming
capabilities. Instead of introducing its own expression syntax, it uses Python
expressions in its preprocessor directives, offering the consistency and
versatility of Python when formulating metaprogramming tasks. It puts strong
emphasis on robustness and on neat integration into developing toolchains.

The project is `hosted on github <https://github.com/aradi/fypp>`_.

`Detailed DOCUMENTATION <http://fypp.readthedocs.org>`_ is available on
`readthedocs.org <http://fypp.readthedocs.org>`_.

Fypp is released under the *BSD 2-clause license*.


Main features
=============

* Definition, evaluation and removal of variables::

    #:if DEBUG > 0
      print *, "Some debug information"
    #:endif

    #:set LOGLEVEL = 2
    print *, "LOGLEVEL: ${LOGLEVEL}$"

    #:del LOGLEVEL

* Macro definitions and macro calls::

    #:def ASSERT(cond)
      #:if DEBUG > 0
        if (.not. ${cond}$) then
          print *, "Assert failed in file ${_FILE_}$, line ${_LINE_}$"
          error stop
        end if
      #:endif
    #:enddef ASSERT

    ! Invoked via direct call (argument needs no quotation)
    @:ASSERT(size(myArray) > 0)

    ! Invoked as Python expression (argument needs quotation)
    $:ASSERT('size(myArray) > 0')

* Conditional output::

    program test
    #:if defined('WITH_MPI')
      use mpi
    #:elif defined('WITH_OPENMP')
      use openmp
    #:else
      use serial
    #:endif

* Iterated output (e.g. for generating Fortran templates)::

    interface myfunc
    #:for dtype in ['real', 'dreal', 'complex', 'dcomplex']
      module procedure myfunc_${dtype}$
    #:endfor
    end interface myfunc

* Inline directives::

    logical, parameter :: hasMpi = #{if defined('MPI')}# .true. #{else}# .false. #{endif}#

* Insertion of arbitrary Python expressions::

    character(*), parameter :: comp_date = "${time.strftime('%Y-%m-%d')}$"

* Inclusion of files during preprocessing::

    #:include "macrodefs.fypp"

* Using Fortran-style continutation lines in preprocessor directives::

    #:if var1 > var2 &
        & or var2 > var4
      print *, "Doing something here"
    #:endif

* Passing (unquoted) multiline string arguments to callables::

    #! Callable needs only string argument
    #:def DEBUG_CODE(code)
      #:if DEBUG > 0
        $:code
      #:endif
    #:enddef DEBUG_CODE

    #! Pass code block as first positional argument
    #:block DEBUG_CODE
      if (size(array) > 100) then
        print *, "DEBUG: spuriously large array"
      end if
    #:endblock DEBUG_CODE

    #! Callable needs also non-string argument types
    #:def REPEAT_CODE(code, repeat)
      #:for ind in range(repeat)
        $:code
      #:endfor
    #:enddef REPEAT_CODE

    #! Pass code block as positional argument and 3 as keyword argument "repeat"
    #:block REPEAT_CODE(repeat=3)
    this will be repeated 3 times
    #:endblock REPEAT_CODE

* Preprocessor comments::

    #! This will not show up in the output
    #! Also the newline characters at the end of the lines will be suppressed

* Suppressing the preprocessor output in selected regions::

    #! Definitions are read, but no output (e.g. newlines) will be produced
    #:mute
    #:include "macrodefs.fypp"
    #:endmute

* Explicit request for stopping the preprocessor::

    #:if DEBUGLEVEL < 0
      #:stop 'Negative debug level not allowed!'
    #:endif

* Easy check for macro parameter sanity::

    #:def mymacro(RANK)
      #! Macro only works for RANK 1 and above
      #:assert RANK > 0
      :
    #:enddef mymacro

* Line numbering directives in output::

    program test
    #:if defined('MPI')
    use mpi
    #:endif
    :

  transformed to ::

    # 1 "test.fypp" 1
    program test
    # 3 "test.fypp"
    use mpi
    # 5 "test.fypp"
    :

  when variable ``MPI`` is defined and Fypp was instructed to generate line
  markers.

* Automatic folding of generated lines exceeding line length limit


Installing
==========

Fypp needs a working Python interpreter. It is compatible with Python 2 (version
2.6 and above) and Python 3 (all versions).

Automatic install
-----------------

Use Pythons command line installer ``pip`` in order to download the stable
release from the `Fypp page on PyPI <http://pypi.python.org/pypi/fypp>`_ and
install it on your system::

  pip install fypp

This installs both, the command line tool ``fypp`` and the Python module
``fypp.py``. Latter you can import if you want to access the functionality of
Fypp directly from within your Python scripts.


Manual install
--------------

For a manual install, you can download the source code of the **stable**
releases from the `Fypp project website
<https://github.com/aradi/fypp/releases>`_.

If you wish to obtain the latest **development** version, clone the projects
repository::

  git clone https://github.com/aradi/fypp.git

and check out the `master` branch.

The command line tool is a single stand-alone script. You can run it directly
from the source folder ::

  FYPP_SOURCE_FOLDER/bin/fypp

or after copying it from the `bin` folder to any location listed in your `PATH`
environment variable, by just issuing ::

  fypp

The python module ``fypp.py`` can be found in ``FYP_SOURCE_FOLDER/src``.


Running
=======

The Fypp command line tool reads a file, preprocesses it and writes it to
another file, so you would typically invoke it like::

  fypp source.fpp source.f90

which would process `source.fpp` and write the result to `source.f90`.  If
input and output files are not specified, information is read from stdin and
written to stdout.

The behavior of Fypp can be influenced with various command line options. A
summary of all command line options can be obtained by::

  fypp -h
