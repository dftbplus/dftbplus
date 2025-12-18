!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_mixer_simple
  use fortuno_serial, only : all_close, suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use dftbp_mixer_simplemixer, only : TSimpleMixerInp, TSimpleMixerReal, TSimpleMixerReal_init
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains

  ! Checks the correctness of the simple mixer's mix1D routine.
  ! The simple mixer works by simply adding chargeDifference weighed by mixParam to the existing
  ! charges.
  $:TEST("simple_mixer_mix1D")
    type(TSimpleMixerInp) :: mixerinp
    type(TSimpleMixerReal) :: mixer

    real(dp), parameter :: atol = 100.0_dp * epsilon(0.0_dp)
    real(dp), parameter :: chargeDiff(3)  = [1.0_dp, 0.0_dp, 2.0_dp]
    real(dp), parameter :: expected(3) = [0.5_dp, 1.0_dp, 3.0_dp]

    real(dp) :: charges(3) 

    ! Prepare mixer
    mixerinp%mixParam = 0.5_dp
    call TSimpleMixerReal_init(mixer, mixerinp)
    call mixer%reset(3)
    
    ! Initial charge configuration
    charges = [0.0_dp, 1.0_dp, 2.0_dp]

    call mixer%mix1D(charges, chargeDiff)
    @:ASSERT(all_close(charges, expected, atol=atol))
  $:END_TEST()

  ! Checks if mix3d properly maps the data to 1D and subsequently invokes the mix1D routine.
  ! Uses the simple mixer, which just add the charge difference weighted by the mix parameter for
  ! simplicity.
  $:TEST("simple_mixer_mix3d")
    type(TSimpleMixerInp) :: mixerinp
    type(TSimpleMixerReal) :: mixer

    real(dp), parameter :: atol = 100.0_dp * epsilon(0.0_dp)
    real(dp) :: charges(3, 3, 3)
    real(dp) :: chargeDiff(3, 3, 3)
    real(dp) :: expected(3, 3, 3)

    ! Since these matrices become quite unwieldy,
    ! we initialise them to 0 and just add a few markers values here and there.
    charges = 0.0_dp 
    chargeDiff = 0.0_dp
    expected = 0.0_dp

    charges(1, 1, 1) = 0.0_dp
    chargeDiff(1, 1, 1) = 1.0_dp
    expected(1, 1, 1) = 0.5_dp

    charges(2, 2, 2) = 1.0_dp
    chargeDiff(2, 2, 2) = 0.0_dp
    expected(2, 2, 2) = 1.0_dp
    
    charges(3, 3, 3) = 10.0_dp
    chargeDiff(3, 3, 3) = 20.0_dp
    expected(3, 3, 3) = 20.0_dp
    
    ! Prepare mixer
    mixerinp%mixParam = 0.5_dp
    call TSimpleMixerReal_init(mixer, mixerinp)
    call mixer%reset(27)

    ! Mix and check
    call mixer%mix3D(charges, chargeDiff)
    @:ASSERT(all_close(reshape(charges, [size(charges)]), reshape(expected, [size(expected)]),&
        & atol=atol))

  $:END_TEST()

  ! Checks if mix6d properly maps the data to 1D and subsequently invokes the mix1D routine.
  $:TEST("simple_mixer_mix6d")
    type(TSimpleMixerInp) :: mixerinp
    type(TSimpleMixerReal) :: mixer

    real(dp), parameter :: atol = 100.0_dp * epsilon(0.0_dp)
    real(dp) :: charges(3, 3, 3, 3, 3, 3)
    real(dp) :: chargeDiff(3, 3, 3, 3, 3, 3)
    real(dp) :: expected(3, 3, 3, 3, 3, 3)
    
    charges = 0.0_dp
    chargeDiff = 0.0_dp
    expected = 0.0_dp

    ! Some markers to catch any obvious mapping errors
    charges(1, 1, 1, 1, 1, 1) = 1.0_dp
    chargeDiff(1, 1, 1, 1, 1, 1) = 1.0_dp
    expected(1, 1, 1, 1, 1 , 1) = 1.5_dp

    charges(2, 2, 2, 2, 2, 2) = 1.0_dp
    chargeDiff(2, 2, 2, 2, 2, 2) = 0.0_dp
    expected(2, 2, 2, 2, 2, 2) = 1.0_dp

    charges(3, 3, 3, 3, 3, 3) = 10.0_dp
    chargeDiff(3, 3, 3, 3, 3, 3) = 20.0_dp
    expected(3, 3, 3, 3, 3, 3) = 20.0_dp

    ! Prepare mixer
    mixerinp%mixParam = 0.5_dp
    call TSimpleMixerReal_init(mixer, mixerinp)
    call mixer%reset(729)

    ! Mix and check
    call mixer%mix6D(charges, chargeDiff)
    @:ASSERT(all_close(reshape(charges, [size(charges)]), reshape(expected, [size(expected)]),&
        & atol=atol))
  $:END_TEST()

  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("simple", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_mixer_simple
