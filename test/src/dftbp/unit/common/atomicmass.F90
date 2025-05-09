!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_common_atomicmass
  use fortuno_serial, only : suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use dftbp_common_atomicmass, only : getAtomicSymbol
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("testSomeElements")
    @:ASSERT(getAtomicSymbol(1.0_dp) == 'H')
    @:ASSERT(getAtomicSymbol(12.0_dp) == 'C')
    @:ASSERT(getAtomicSymbol(16.0_dp) == 'O')
  $:END_TEST()


  $:TEST("testTooLargeMass")
    @:ASSERT(getAtomicSymbol(500.0_dp) == '??')
  $:END_TEST()


  $:TEST("testNegativeMass")
    @:ASSERT(getAtomicSymbol(-1.0_dp) == '??')
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("atomicmass", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_common_atomicmass
