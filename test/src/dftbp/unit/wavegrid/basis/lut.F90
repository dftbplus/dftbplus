!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_wavegrid_basis_lut
  use fortuno_serial, only : is_close, suite => serial_suite_item, test_list
  use dftbp_wavegrid_basis, only : TRadialTableOrbital, TSlaterOrbital, TRadialTableOrbital_initFromArray, &
       TRadialTableOrbital_initFromOrbital, TSlaterOrbital_init
  use dftbp_common_accuracy, only : dp
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none
  
  private
  public :: tests

  !> Allow 0.001% relative error
  real(dp), parameter :: rtol = 1.0e-5_dp


contains
  !> Check if initialisation from existing lookup table and subsequent 
  !! Access works as expected.
  $:TEST("TRadialTableOrbital_initFromArray")
    type(TRadialTableOrbital) :: sto
    real(dp) :: r, val, expected
    integer :: i
    integer, parameter :: angMom = 0 ! Only relevant for realTessY
    real(dp), parameter :: gridDist = 1.0_dp
    ! f(x)=x^2
    real(dp), parameter :: gridValue(5) = [0.0_dp, 1.0_dp, 4.0_dp, 9.0_dp, 16.0_dp]
    ! Lookup table access is simple linear interpolation 
    real(dp), parameter :: checkedPairs(2,7) = reshape([&
        0.0_dp, 0.0_dp, &
        0.5_dp, 0.5_dp, &
        1.0_dp, 1.0_dp, &
        1.5_dp, 2.5_dp, &
        3.1_dp, 0.9_dp*9.0_dp + 0.1_dp*16.0_dp, &
        4.0_dp, 0.0_dp, &  ! no i+1 for interpolation
        4.01_dp, 0.0_dp &  ! beyond cutoff
      ], [2,7])


    call TRadialTableOrbital_initFromArray(sto, gridValue, gridDist, angMom)
    @:ASSERT(sto%angMom == angMom)

    do i = 1, size(checkedPairs, 2)
      r = checkedPairs(1, i)
      expected = checkedPairs(2, i)
      val = sto%getRadial(r)
      @:CHECK(is_close(val, expected, rtol=rtol))
    end do 
  $:END_TEST()


  !> Check resampling of existing analytical orbital to lookup table.
  $:TEST("TRadialTableOrbital_initFromOrbital")
    type(TSlaterOrbital) :: sto_1s
    type(TRadialTableOrbital) :: lut
    real(dp), parameter :: aa(1,1) = reshape([2.0_dp], [1,1])
    real(dp), parameter :: alpha(1) = [1.0_dp]
    real(dp) :: resolution, r, valSto, valLut

    call TSlaterOrbital_init(sto_1s, aa, alpha, angMom=0, cutoff=15.0_dp)

    resolution = 0.01_dp
    call TRadialTableOrbital_initFromOrbital(lut, sto_1s, resolution)

    ! Check if an interpolated point is close to the original function
    r = 3.5_dp * resolution
    valSto = sto_1s%getRadial(r)
    valLut = lut%getRadial(r)
    @:CHECK(is_close(valSto, valLut, rtol=1.0e-4_dp))
  $:END_TEST()


  !> Register test cases with Fortuno
  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("lut", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_wavegrid_basis_lut
