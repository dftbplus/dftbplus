!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_mixer_diis
  use fortuno_serial, only : all_close, suite => serial_suite_item, test_list
  use dftbp_common_accuracy, only : dp
  use dftbp_mixer_diismixer, only : TDiisMixerInp, TDiisMixerReal, TDiisMixerReal_init
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains

  ! Basic Regression test to notify about potentially unintended changes to the DIIS mixer behavior.
  ! Charges and Differences obtained from a simple H3 calculation.
  $:TEST("diis_mixer_H3")
    type(TDiisMixerInp) :: inp
    type(TDiisMixerReal) :: mixer
    
    real(dp), parameter :: atol = 100.0_dp * epsilon(0.0_dp)
    real(dp), parameter :: expectedCharges(3, 6) = reshape([ &
        1.0000000000000000_dp, 1.0000000000000000_dp, 1.0000000000000000_dp, &
        0.76814995249643614_dp, 1.3486921095027322_dp, 0.88315793800083164_dp, &
        0.84162976921716104_dp, 1.2336229862547019_dp, 0.92474724452813706_dp, &
        0.82851213597245810_dp, 1.2547909429916608_dp, 0.91669692103588141_dp, &
        0.82855095868412298_dp, 1.2547154075907785_dp, 0.91673363372509797_dp, &
        0.82855125884896574_dp, 1.2547148326356727_dp, 0.91673390851536196_dp], [3,6]) 
    real(dp), parameter :: difference(3, 5) = reshape([ &
        -0.19320837291963655_dp, 0.29057675791894355_dp, -9.7368384999307001E-002_dp, &
        6.7711031765136043E-002_dp, -0.10672264480614602_dp, 3.9011613041009863E-002_dp, &
        -1.4661076723682331E-002_dp, 2.3994941612978993E-002_dp, -9.3338648892962173E-003_dp, &
        4.3509104117211983E-005_dp, -8.6630676363164127E-005_dp, 4.3121572245508055E-005_dp, &
        3.3405783306239556E-007_dp, -6.5440257346338626E-007_dp, 3.2034474140019142E-007_dp], [3,5])

    real(dp) :: qInpResult(3) 
    integer :: ii

    ! Prepare mixer
    inp%iGenerations = 3
    inp%initMixParam = 0.2_dp
    inp%tFromStart = .true.
    inp%alpha = 0.0_dp
    call TDiisMixerReal_init(mixer, inp)
    call mixer%reset(3)

    ! Assert that each pair of (charges, chargeDifference) results in the expected charges.
    do ii = 1, 5
        qInpResult = expectedCharges(:,ii)
        call mixer%mix(qInpResult, difference(:,ii))
        @:ASSERT(all_close(qInpResult, expectedCharges(:,ii+1), atol=atol))
    end do
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("diis", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_mixer_diis
