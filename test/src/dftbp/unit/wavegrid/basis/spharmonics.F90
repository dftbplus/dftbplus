!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_wavegrid_basis_spharmonics
  use fortuno_serial, only : is_close, suite => serial_suite_item, test_list
  use test_wavegrid_basis_spharmonics_ref, only : spotCheck, spotChecks_l0, spotChecks_l1, &
    & spotChecks_l2, spotChecks_l3, spotChecks_l4, atol
  use dftbp_wavegrid_basis, only : realTessY
  use dftbp_common_accuracy, only : dp
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none
  
  private
  public :: tests

contains

  subroutine checkBatch(spotChecks)
    type(spotCheck), intent(in) :: spotChecks(:)
    real(dp) :: ylm, inv_r
    type(spotCheck) :: c
    integer :: ii

    do ii = 1, size(spotChecks)
      c = spotChecks(ii)
      inv_r = 1.0_dp / sqrt(c%x**2 + c%y**2 + c%z**2)
      ylm = realTessY(c%l, c%m, [c%x, c%y, c%z], inv_r)
      @:CHECK(is_close(ylm, c%expected, atol=atol))
    end do 
  end subroutine checkBatch

  $:TEST("spharmonicsl0")
    call checkBatch(spotChecks_l0)
  $:END_TEST()

  $:TEST("spharmonicsl1")
    call checkBatch(spotChecks_l1)
  $:END_TEST()

  $:TEST("spharmonicsl2")
    call checkBatch(spotChecks_l2)
  $:END_TEST()

  $:TEST("spharmonicsl3")
    call checkBatch(spotChecks_l3)
  $:END_TEST()

  $:TEST("spharmonicsl4")
    call checkBatch(spotChecks_l4)
  $:END_TEST()


  !> Register test cases with Fortuno
  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("spharmonics", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_wavegrid_basis_spharmonics
