!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_serial.fypp"

module test_elecsolvers_partialdiag
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_etemp, only : fillingTypes
  use dftbp_elecsolvers_partialdiag, only : getEmptyStateWindow, getNOccupied, TPartialDiag,&
      & TPartialDiag_init
  use fortuno_serial, only : suite => serial_suite_item, test_list
  $:FORTUNO_SERIAL_IMPORTS()
  implicit none

  private
  public :: tests

contains


  $:TEST("occupied_states_closed_shell")
    ! two electrons per state
    @:ASSERT(getNOccupied([8.0_dp], 1) == 4)
    ! an odd number of electrons needs a further, partly filled, state
    @:ASSERT(getNOccupied([9.0_dp], 1) == 5)
  $:END_TEST()


  $:TEST("occupied_states_collinear_spin")
    ! the channels are diagonalised separately, so the larger one sets the count
    @:ASSERT(getNOccupied([5.0_dp, 3.0_dp], 2) == 5)
    @:ASSERT(getNOccupied([3.0_dp, 5.0_dp], 2) == 5)
  $:END_TEST()


  $:TEST("occupied_states_pauli")
    ! one electron per state of the combined problem
    @:ASSERT(getNOccupied([8.0_dp], 4) == 8)
  $:END_TEST()


  $:TEST("empty_state_window")
    real(dp), parameter :: kT = 0.01_dp
    real(dp), parameter :: tol = 1.0e-12_dp

    ! the Fermi-Dirac tail decays more slowly than the Gaussian based ones
    @:ASSERT(abs(getEmptyStateWindow(fillingTypes%Fermi, kT) - 24.0_dp * kT) < tol)
    @:ASSERT(abs(getEmptyStateWindow(fillingTypes%Gaussian, kT) - 5.0_dp * kT) < tol)
    @:ASSERT(abs(getEmptyStateWindow(fillingTypes%Methfessel, kT) - 5.0_dp * kT) < tol)
    @:ASSERT(getEmptyStateWindow(fillingTypes%Fermi, 0.0_dp) < tol)
  $:END_TEST()


  $:TEST("inactive_by_default")
    type(TPartialDiag) :: partialDiag

    ! a negative number of empty states means the whole spectrum
    call TPartialDiag_init(partialDiag, -1, 100, [8.0_dp], 1)
    @:ASSERT(.not. partialDiag%isActive)
    @:ASSERT(partialDiag%getNState() == 100)
    @:ASSERT(partialDiag%isSufficient)
  $:END_TEST()


  $:TEST("requested_states")
    type(TPartialDiag) :: partialDiag

    call TPartialDiag_init(partialDiag, 10, 100, [8.0_dp], 1)
    @:ASSERT(partialDiag%isActive)
    @:ASSERT(partialDiag%getNState() == 14)

    ! more empty states than the basis holds is not possible
    call TPartialDiag_init(partialDiag, 1000, 100, [8.0_dp], 1)
    @:ASSERT(partialDiag%getNState() == 100)
  $:END_TEST()


  $:TEST("growth_doubles_and_saturates")
    type(TPartialDiag) :: partialDiag

    call TPartialDiag_init(partialDiag, 10, 100, [8.0_dp], 1)
    call partialDiag%grow()
    @:ASSERT(partialDiag%getNState() == 24)
    call partialDiag%grow()
    @:ASSERT(partialDiag%getNState() == 44)

    ! growth stops at a full diagonalisation
    call partialDiag%grow()
    call partialDiag%grow()
    call partialDiag%grow()
    @:ASSERT(partialDiag%getNState() == 100)

    ! a single empty state doubles like any other count
    call TPartialDiag_init(partialDiag, 1, 100, [8.0_dp], 1)
    @:ASSERT(partialDiag%getNState() == 5)
    call partialDiag%grow()
    @:ASSERT(partialDiag%getNState() == 6)
  $:END_TEST()


  $:TEST("sufficiency_of_calculated_states")
    type(TPartialDiag) :: partialDiag
    real(dp) :: eigvals(6, 1, 1)
    real(dp), parameter :: kT = 0.01_dp
    real(dp), parameter :: Ef(1) = [0.0_dp]

    call TPartialDiag_init(partialDiag, 2, 100, [8.0_dp], 1)

    ! highest calculated state well above the window, so anything left out is empty
    eigvals(:, 1, 1) = [-1.0_dp, -0.9_dp, -0.8_dp, -0.7_dp, -0.6_dp, 1.0_dp]
    call partialDiag%check(eigvals, Ef, kT, fillingTypes%Fermi)
    @:ASSERT(partialDiag%isSufficient)

    ! highest calculated state inside the window, so states above it may still be occupied
    eigvals(6, 1, 1) = 0.5_dp * 24.0_dp * kT
    call partialDiag%check(eigvals, Ef, kT, fillingTypes%Fermi)
    @:ASSERT(.not. partialDiag%isSufficient)

    ! the same spectrum is sufficient for a Gaussian filling, whose tail is shorter
    eigvals(6, 1, 1) = 6.0_dp * kT
    call partialDiag%check(eigvals, Ef, kT, fillingTypes%Gaussian)
    @:ASSERT(partialDiag%isSufficient)
    call partialDiag%check(eigvals, Ef, kT, fillingTypes%Fermi)
    @:ASSERT(.not. partialDiag%isSufficient)
  $:END_TEST()


  $:TEST("sufficiency_across_k_points_and_spins")
    type(TPartialDiag) :: partialDiag
    real(dp) :: eigvals(3, 2, 2)
    real(dp), parameter :: kT = 0.01_dp
    real(dp), parameter :: Ef(2) = [0.0_dp, 0.0_dp]

    call TPartialDiag_init(partialDiag, 1, 100, [4.0_dp, 4.0_dp], 2)

    ! every k-point and spin channel has its own spectrum, all of them need enough empty states
    eigvals(:,:,:) = 1.0_dp
    call partialDiag%check(eigvals, Ef, kT, fillingTypes%Fermi)
    @:ASSERT(partialDiag%isSufficient)

    eigvals(3, 2, 1) = 0.01_dp * kT
    call partialDiag%check(eigvals, Ef, kT, fillingTypes%Fermi)
    @:ASSERT(.not. partialDiag%isSufficient)
  $:END_TEST()


  $:TEST("full_spectrum_is_always_sufficient")
    type(TPartialDiag) :: partialDiag
    real(dp) :: eigvals(6, 1, 1)
    real(dp), parameter :: Ef(1) = [1.0e3_dp]

    ! nothing was left out, so nothing can be missing, whatever the eigenvalues are
    call TPartialDiag_init(partialDiag, 2, 6, [8.0_dp], 1)
    eigvals(:, 1, 1) = 0.0_dp
    call partialDiag%check(eigvals, Ef, 0.01_dp, fillingTypes%Fermi)
    @:ASSERT(partialDiag%isSufficient)
  $:END_TEST()


  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("partialdiag", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_elecsolvers_partialdiag
