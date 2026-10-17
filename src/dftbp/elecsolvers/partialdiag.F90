!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Bookkeeping for partial diagonalisation, where the eigensolver is only asked for the occupied
!! states together with a limited number of empty states.
!!
!! The number of empty states is chosen by the user as a starting value and then increased during
!! the self-consistent cycle whenever the states obtained are found to be insufficient. The test for
!! sufficiency is that the highest calculated state lies far enough above the chemical potential for
!! its occupation to be numerically negligible, the distance being a multiple of the broadening
!! width of the filling function in use.
module dftbp_elecsolvers_partialdiag
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_etemp, only : fillingTypes
  implicit none

  private
  public :: TPartialDiag, TPartialDiag_init
  public :: getEmptyStateWindow, getNOccupied


  !> Occupation small enough to be neglected when deciding whether a state is empty. It has to be
  !! well below the accuracy to which the charges are converged, as the neglected states would
  !! otherwise leave a systematic error behind which the self-consistent cycle cannot remove.
  real(dp), parameter :: negligibleFilling = 1.0e-10_dp

  !> Multiple of the broadening width above the chemical potential at which the Fermi-Dirac
  !! occupation 1 / (1 + exp(x)) falls below negligibleFilling, which is ln(1 / f - 1) = 23.0,
  !! rounded up.
  real(dp), parameter :: nWidthsFermi = 24.0_dp

  !> Same for the Gaussian and Methfessel-Paxton schemes. Their occupations are integrals of the
  !! broadening function of the level, so the point where the broadening function exp(-x**2) itself
  !! falls below negligibleFilling is sufficient, which is sqrt(ln(1 / f)) = 4.8, rounded up. The
  !! Gaussian is the slowest decaying member of the family, so this covers the higher orders too.
  real(dp), parameter :: nWidthsGaussian = 5.0_dp


  !> State counting for partial diagonalisation
  type :: TPartialDiag

    !> Whether only a part of the spectrum is being calculated
    logical :: isActive = .false.

    !> Number of empty states currently being requested
    integer :: nEmpty = 0

    !> Number of states needed to hold the electrons of the system
    integer :: nOccupied = 0

    !> Number of states in the basis, i.e. the size of a full diagonalisation
    integer :: nOrb = 0

    !> Whether the last diagonalisation returned at least one completely empty state
    logical :: isSufficient = .true.

  contains

    procedure :: getNState => TPartialDiag_getNState
    procedure :: check => TPartialDiag_check
    procedure :: grow => TPartialDiag_grow

  end type TPartialDiag


contains


  !> Initialises the state counting.
  subroutine TPartialDiag_init(this, nEmpty, nOrb, nEl, nSpin)

    !> Instance
    type(TPartialDiag), intent(out) :: this

    !> Number of empty states to start from, negative for a full diagonalisation
    integer, intent(in) :: nEmpty

    !> Number of states in the basis
    integer, intent(in) :: nOrb

    !> Number of electrons, per spin channel where these are separate
    real(dp), intent(in) :: nEl(:)

    !> Number of spin blocks in the Hamiltonian
    integer, intent(in) :: nSpin

    this%nOrb = nOrb
    this%nOccupied = min(getNOccupied(nEl, nSpin), nOrb)
    this%isActive = (nEmpty >= 0)
    if (this%isActive) then
      this%nEmpty = min(nEmpty, nOrb - this%nOccupied)
    else
      this%nEmpty = nOrb - this%nOccupied
    end if
    this%isSufficient = .true.

  end subroutine TPartialDiag_init


  !> Number of states which the eigensolver should return.
  pure function TPartialDiag_getNState(this) result(nState)

    !> Instance
    class(TPartialDiag), intent(in) :: this

    !> Number of states to calculate
    integer :: nState

    nState = min(this%nOccupied + this%nEmpty, this%nOrb)

  end function TPartialDiag_getNState


  !> Tests whether the calculated states extend far enough above the chemical potential that the
  !! states which were not calculated are certainly empty.
  !!
  !! The test is applied to the highest calculated state of every k-point and spin channel, as the
  !! spectra of these are independent of each other.
  subroutine TPartialDiag_check(this, eigvals, Ef, tempElec, iDistribFn)

    !> Instance
    class(TPartialDiag), intent(inout) :: this

    !> Calculated eigenvalues (state, k-point, spin), holding the calculated states only
    real(dp), intent(in) :: eigvals(:,:,:)

    !> Chemical potential(s)
    real(dp), intent(in) :: Ef(:)

    !> Electronic broadening width, in energy units
    real(dp), intent(in) :: tempElec

    !> Filling function in use
    integer, intent(in) :: iDistribFn

    integer :: nState

    if (.not. this%isActive) then
      this%isSufficient = .true.
      return
    end if

    nState = size(eigvals, dim=1)
    if (nState >= this%nOrb) then
      ! nothing was left out, so nothing can be missing
      this%isSufficient = .true.
      return
    end if

    this%isSufficient = minval(eigvals(nState,:,:))&
        & > maxval(Ef) + getEmptyStateWindow(iDistribFn, tempElec)

  end subroutine TPartialDiag_check


  !> Increases the number of empty states, doubling it until the full basis is reached.
  subroutine TPartialDiag_grow(this)

    !> Instance
    class(TPartialDiag), intent(inout) :: this

    if (.not. this%isActive) then
      return
    end if

    this%nEmpty = min(max(2 * this%nEmpty, 1), this%nOrb - this%nOccupied)

  end subroutine TPartialDiag_grow


  !> Distance above the chemical potential beyond which the occupation of a state is numerically
  !! negligible for the given filling function.
  pure function getEmptyStateWindow(iDistribFn, tempElec) result(window)

    !> Filling function in use
    integer, intent(in) :: iDistribFn

    !> Electronic broadening width, in energy units
    real(dp), intent(in) :: tempElec

    !> Resulting energy window
    real(dp) :: window

    if (iDistribFn == fillingTypes%Fermi) then
      window = nWidthsFermi * tempElec
    else
      window = nWidthsGaussian * tempElec
    end if

  end function getEmptyStateWindow


  !> Number of states required to hold the electrons of the system.
  !!
  !! For collinear spin the two channels are diagonalised separately, so the larger of the two
  !! channels determines the count. Where the Fermi level is shared between the channels, electrons
  !! can move between them during the self-consistent cycle, which the adaptive increase of the
  !! number of empty states takes care of.
  pure function getNOccupied(nEl, nSpin) result(nOccupied)

    !> Number of electrons, per spin channel where these are separate
    real(dp), intent(in) :: nEl(:)

    !> Number of spin blocks in the Hamiltonian
    integer, intent(in) :: nSpin

    !> Number of states needed for these electrons
    integer :: nOccupied

    select case (nSpin)
    case (1)
      ! two electrons per state
      nOccupied = ceiling(0.5_dp * sum(nEl))
    case (2)
      ! one electron per state, the channels being diagonalised separately
      nOccupied = ceiling(maxval(nEl))
    case default
      ! one electron per state of the combined Pauli problem
      nOccupied = ceiling(sum(nEl))
    end select

    nOccupied = max(nOccupied, 1)

  end function getNOccupied

end module dftbp_elecsolvers_partialdiag
