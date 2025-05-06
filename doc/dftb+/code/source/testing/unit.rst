*********************
Unit testing in DFTB+
*********************

The 24.1 release candidate, as of commit `62ad90f85` now uses the
`Fortuno <https://fortuno.readthedocs.io/en/latest/>`_ framework for
unit testing. The tests are located in a subdirectory of the `test/`
folder, with a directory structure mirroring the main source tree in
`src/dftbp/` and tests included via cmake::

  test/src/dftbp/unit
  ├── test/src/dftbp/unit/CMakeLists.txt
  ├── test/src/dftbp/unit/common
  ├── test/src/dftbp/unit/dftb
  ├── test/src/dftbp/unit/include
  ├── test/src/dftbp/unit/io
  ├── test/src/dftbp/unit/math
  └── test/src/dftbp/unit/type

Earlier code versions used a different system, with unit test cases
located in `test/src/dftbp/unittests`, so if you are developing code
based on earlier releases and have added unit tests, you will need to
adapt them to the new structure.

The `fypp <https://fypp.readthedocs.io/en/stable/fypp.html>`_ enhanced
version of Fortuno is used, so test cases should include the macros
and use `@:ASSERT()` to test statements that are true about results. A
minimal test case can look something like

.. code-block:: Fortran
   :linenos:

   #:include "fortuno_serial.fypp"

   module test_common_accuracy ! name corresponds to test_ then path to code being tested
     use dftbp_common_accuracy, only : cp, dp ! import what is need for test(s)
     use fortuno_serial, only : suite => serial_suite_item, test_list
     $:FORTUNO_SERIAL_IMPORTS()
     implicit none

     private
     public :: tests

   contains

      $:TEST("trivial")

        integer, parameter :: ii = 0

        ! Actual test, this does something very trivial:
        @:ASSERT(ii == 0)

      $:END_TEST()


      $:TEST("types")

        ! Actual test, this does something very trivial with imported variables:
        @:ASSERT(cp == dp)

      $:END_TEST()


      ! boiler-plate code to run the tests

      function tests()
        type(test_list) :: tests

          tests = test_list([&
              suite("accuracy", test_list([&
                  $:TEST_ITEMS()
              ]))&
          ])
          $:STOP_ON_MISSING_TEST_ITEMS()

      end function tests

   end module test_common_accuracy


See existing test examples in the repository for more realistic cases.
