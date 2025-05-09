!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_mixer_anderson
  use fortuno_serial, only : all_close, suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use dftbp_mixer_andersonmixer, only : TAndersonMixerInp, TAndersonMixerReal,&
      & TAndersonMixerReal_init
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains
    
  ! Basic Regression test to notify about any changes to Anderson mixer behavior.
  ! Charges and Differences obtained from a simple H2O calculation.
  $:TEST("anderson_mixer_H2O")
    type(TAndersonMixerInp) :: inp
    type(TAndersonMixerReal) :: mixer
    
    real(dp), parameter :: atol = 100.0_dp * epsilon(0.0_dp)
    real(dp), parameter :: expectedCharges(3, 5) = reshape([ &
        6.0000000000000018_dp, 1.0000000000000002_dp,  1.0000000000000002_dp,  &
        6.0088081626601531_dp, 0.99559591866992481_dp, 0.99559591866992481_dp, &
        6.4903876826943883_dp, 0.75480615865281897_dp, 0.75480615865281908_dp, &
        6.4622134796623394_dp, 0.76889326016884263_dp, 0.76889326016884274_dp, &
        6.4635958158228206_dp, 0.76820209208860180_dp, 0.76820209208860191_dp], [3,5]) 
    real(dp), parameter :: difference(3, 4) = reshape([ &
       0.88081626601508667_dp, -0.44040813300754278_dp, -0.44040813300754278_dp, &
       0.86499689478346919_dp, -0.43249844739173582_dp, -0.43249844739173560_dp, &
       -5.3722892233180275E-002_dp, 2.6861446116577481E-002_dp, 2.6861446116578147E-002_dp, &
       2.6725939528224885E-003_dp, -1.3362969764226795E-003_dp, -1.3362969764224575E-003_dp], [3,4])

    real(dp) :: qInpResult(3)
    integer :: ii

    ! Prepare mixer
    inp%iGenerations = 4
    inp%mixParam = 0.05_dp
    inp%initMixParam = 0.01_dp
    call TAndersonMixerReal_init(mixer, inp)
    call mixer%reset(3)
   
    ! Load initial charges
    qInpResult = expectedCharges(:,1)

    ! Assert that each pair of (charges, chargeDifference) results in the expected charges.
    do ii = 1, 4
        call mixer%mix(qInpResult, difference(:,ii))
        @:ASSERT(all_close(qInpResult, expectedCharges(:,ii+1), atol=atol))
    end do
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("anderson", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_mixer_anderson
