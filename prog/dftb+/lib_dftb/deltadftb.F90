!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines for (time independent excited) TI-DFTB
module dftbp_deltadftb
  use dftbp_accuracy, only : dp
  use dftbp_energytypes, only : TEnergies
  implicit none

  private
  public :: ZieglerSum, determinants, initTIDFTB

  !> Name space for determinants
  type :: TDeterminantsEnum

    integer :: ground = 0
    integer :: triplet = 1
    integer :: mixed = 2

  end type TDeterminantsEnum

  !> Actual values for elecSolverTypes.
  type(TDeterminantsEnum), parameter :: determinants = TDeterminantsEnum()

contains

  !> Spin Purifies Non-Aufbau excited state energy and forces
  subroutine ZieglerSum(energy, derivs, tripletderivs, mixedderivs, tGeoOpt)

    !> energy components
    type(TEnergies), intent(inout) :: energy

    !> derivative components
    real(dp), intent(inout), allocatable, optional:: derivs(:,:)

    !> Triplet state derivatives
    real(dp), intent(inout), allocatable, optional :: tripletderivs(:,:)

    !> Spin contaminated derivatives
    real(dp), intent(inout), allocatable, optional :: mixedderivs(:,:)

    !> Is geometry being optimised
    logical, intent(in) :: tGeoOpt


    ! Ziegler sum rule: E_S1 = 2E_mix - E_triplet
    energy%Etotal = 2.0_dp * energy%Emixed - energy%Etriplet
    energy%EMermin = 2.0_dp * energy%EmixMermin - energy%EtripMermin
    energy%Ezero = 2.0_dp * energy%EmixZero - energy%EtripZero
    energy%EGibbs = 2.0_dp * energy%EmixGibbs - energy%EtripGibbs
    energy%EForceRelated = 2.0_dp * energy%EmixForceRelated - energy%EtripForceRelated
    ! dE_S1 = 2dE_mix - dE_triplet
    if (tGeoOpt) then
      derivs(:,:) = 2.0_dp * mixedderivs - tripletderivs
    endif

  end subroutine ZieglerSum


  !> Initialised Time-independent excited state DFTB (TI-DFTB) conditions for determinant
  subroutine initTIDFTB(tNonAufbau, tSpinPurify, tGroundGuess, firstDeterminant, lastDeterminant)

    !> Is this a non-aufbau filled TI calculation
    logical, intent(in) :: tNonAufbau

    !> Is this a spin purified calculation?
    logical, intent(in) :: tSpinPurify

    !> Should there be a ground state intial guess before Non-Aufbau calc?
    logical, intent(in) :: tGroundGuess

    !> First determinant in space to calculate
    integer, intent(out) :: firstDeterminant

    !> Last state being calculated
    integer, intent(out) :: lastDeterminant

    if (.not.tNonAufbau) then
      firstDeterminant = 1
      lastDeterminant = 1
      return
    end if

    if (tGroundGuess) then
      firstDeterminant = determinants%ground
    else
      firstDeterminant = determinants%triplet
    end if
    if (tSpinPurify) then
      lastDeterminant = determinants%mixed
    else
      lastDeterminant = determinants%triplet
    end if

  end subroutine initTIDFTB


end module dftbp_deltadftb
