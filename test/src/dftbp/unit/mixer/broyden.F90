!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_mixer_broyden
  use dftbp_common_accuracy, only : dp
  use dftbp_mixer_broydenmixer, only : TBroydenMixerInp, TBroydenMixerReal, TBroydenMixerReal_init
  use fortuno_serial, only : suite => serial_suite_item, test_list, all_close
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains

  ! Basic regression test to notify about any changes to Broyden mixer behavior.
  ! Charges and Differences obtained from a simple water molecule calculation.
  $:TEST("broyden_mixer_H20")
    type(TBroydenMixerInp) :: inp
    type(TBroydenMixerReal) :: mixer
    
    real(dp), parameter :: atol = 100.0_dp * epsilon(0.0_dp)
    real(dp), parameter :: expectedCharges(3, 5) = reshape([ &
        6.0000000000000018_dp, 1.0000000000000002_dp, 1.0000000000000002_dp, &  
        6.1761632532030193_dp, 0.91191837339849169_dp, 0.91191837339849169_dp, & 
        6.4798000093100434_dp, 0.76009999534497941_dp, 0.76009999534497941_dp, & 
        6.4625859389302871_dp, 0.76870703053485734_dp, 0.76870703053485745_dp, & 
        6.4635561795150211_dp, 0.76822191024249031_dp, 0.76822191024249042_dp], [3,5])
    real(dp), parameter :: difference(3, 4) = reshape([ &
        0.88081626601508667_dp, -0.44040813300754278_dp, -0.44040813300754278_dp, &
        0.55742893386488124_dp, -0.27871446693244140_dp, -0.27871446693244140_dp, &
        -3.2497351757903914E-002_dp, 1.6248675878950958E-002_dp, 1.6248675878951291E-002_dp, &
        1.9288771697691942E-003_dp, -9.6443858488526324E-004_dp, -9.6443858488515222E-004_dp], [3,4])

    real(dp) :: qInpResult(3)
    integer :: ii

    ! Prepare mixer
    inp%mixParam = 0.2_dp
    inp%omega0 = 0.01_dp
    inp%minWeight = 1.0_dp
    inp%maxWeight = 100000.0_dp
    inp%maxSccIter = 100
    inp%weightFac = 0.01_dp
    call TBroydenMixerReal_init(mixer, inp)
    call mixer%reset(3)

    ! Set initial charges
    qInpResult = expectedCharges(:,1)

    ! Iteratively assert that each pair of (charges, chargeDifference) results in the expected charges.
    do ii = 1, 4
        call mixer%mix(qInpResult, difference(:,ii))
        @:ASSERT(all_close(qInpResult, expectedCharges(:,ii+1), atol=atol))
    end do
  $:END_TEST()

  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("broyden", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_mixer_broyden
