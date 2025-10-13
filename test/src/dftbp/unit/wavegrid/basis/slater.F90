!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_wavegrid_basis_slater
  use fortuno_serial, only : is_close, suite => serial_suite_item, test_list
  use dftbp_wavegrid_basis, only : TSlaterOrbital
  use dftbp_common_accuracy, only : dp
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none
  
  private
  public :: tests

  !> Allow 0.001% relative error
  real(dp), parameter :: rtol = 1.0e-5_dp


contains


  $:TEST("TSlaterOrbital_getRadial")
    type(TSlaterOrbital) :: sto_1s, sto_2p
    real(dp) :: cutoff, val, expected
    real(dp), parameter :: aa(1,1) = reshape([2.0_dp], [1,1])
    real(dp), parameter :: alpha(1) = [1.0_dp]

    cutoff = 20.0_dp

    call sto_1s%init(aa=aa, alpha=alpha, angMom=0, cutoff=cutoff)
    val = sto_1s%getRadial(1.0_dp)
    expected = 2.0_dp * exp(-1.0_dp)
    @:CHECK(is_close(val, expected, rtol=rtol))

    call sto_2p%init(aa=reshape([1.0_dp], [1,1]), alpha=alpha, angMom=1, cutoff=cutoff)
    val = sto_2p%getRadial(1.0_dp)
    expected = 1.0_dp * exp(-1.0_dp)
    @:CHECK(is_close(val, expected, rtol=rtol))

    val = sto_2p%getRadial(0.0_dp)
    expected = 0.0_dp
    @:CHECK(is_close(val, expected, atol=1e-15_dp))
  $:END_TEST()


  !> Register test cases with Fortuno
  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("slater", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_wavegrid_basis_slater
