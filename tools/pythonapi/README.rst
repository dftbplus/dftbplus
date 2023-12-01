PYTHONAPI: A ctypes based Python interface for DFTB+
****************************************************

This package provides a Python Interface for DFTB+ to improve
the ease of use and expand its applications and functionality.
It currently contains the following Python class:

DftbPlus
  Interface module for the communication between DFTB+ and
  Python (via the foreign function C-library ctypes). Provides
  methods for initializing and configuring a DFTB+ calculator.


Compiling DFTB+
===============

In order to be able to use the Python interface, DFTB+ has to be
compiled as a shared library with API support enabled. To instruct
cmake that a dynamically linked shared library should be created
containing DFTB+. This can be done by setting the WITH_API, WITH_PYTHON,
BUILD_SHARED_LIBS and ENABLE_DYNAMIC_LOADING flags to be TRUE in the
configuration file config.cmake, before starting the configuration
and compilation process. Alternatively, if you do not want to
modify files, a construct like the following is a convenient way
to specify these flags on the command line while configuring with
CMake:

cmake -DENABLE_DYNAMIC_LOADING=1 -DBUILD_SHARED_LIBS=1 -DWITH_API=1 -DWITH_PYTHON=1 ..


Testing pythonapi
=================

In the _build/ directory, running

ctest -R pyapi_*

will validate the compiled library and the source of the pythonapi by
executing regression tests.

For developers
--------------

To perform pylint static checking from the top level directory of the
DFTB+ project, use

pylint3 --rcfile utils/srccheck/pylint/pylintrc-3.ini tools/pythonapi/src/*


Installation
============

Please note, that this package has been tested for **Python 3.X**
support. It additionally needs Numerical Python (the numpy module).

You can install the script package via the standard 'pip'
mechanism::

  python -m pip install .

with an appropriate level of permission.
