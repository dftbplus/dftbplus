!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines for (time independent excited) TI-DFTB
module dftbp_dftbdeterminants
  use dftbp_accuracy, only : dp
  use dftbp_energytypes, only : TEnergies
  use dftbp_message, only : error
  use dftbp_etemp, only : Efilling
  implicit none

  private
  public :: determinants, TDftbDeterminants_init, TDftbDeterminants

  !> Name space for determinants
  type :: TDeterminantsEnum

    !> Conventional DFTB ground state
    integer :: ground = 0

    !> Triplet configuration
    integer :: triplet = 1

    !> Hole under occupied levels, spin contaminated
    integer :: mixed = 2

  end type TDeterminantsEnum

  !> Actual values for elecSolverTypes.
  type(TDeterminantsEnum), parameter :: determinants = TDeterminantsEnum()

  !> Control type for Delta DFTB / TI-DFTB
  type TDftbDeterminants

    !> Is this a non-Aufbau filling
    logical :: isNonAufbau

    !> Should the non-Aufbau fillings be spin purified
    logical :: isSpinPurify

    !> Should a ground state case be calculated as well
    logical :: isGroundGuess

    !> Current determinant being solved
    integer :: iDeterminant

    !> Number of electrons in each spin channel
    real(dp), allocatable :: nEl(:)

    !> Has the calculation finished and results are now ready to use
    logical :: isFinished

  contains

    procedure :: postProcessDets
    procedure :: nDeterminant
    procedure :: whichDeterminant
    procedure :: detFilling

  end type TDftbDeterminants


contains

  !> Counts number of determiants to be evaluated
  pure function nDeterminant(this) result(nDet)

    !> Instance
    class(TDftbDeterminants), intent(in) :: this

    integer :: nDet

    nDet = 1
    if (.not.this%isNonAufbau) then
      return
    end if

    if (this%isSpinPurify) then
      nDet = 2
    end if

    if (this%isGroundGuess) then
      nDet = nDet + 1
    end if

  end function nDeterminant


  !> Converts determinant number into what type it should be
  function whichDeterminant(this, iDet) result(det)

    !> Instance
    class(TDftbDeterminants), intent(inout) :: this

    integer, intent(in) :: iDet

    integer :: det

    this%isFinished = .false.

    if (.not.this%isNonAufbau) then
      if (iDet /= 1) then
        call error("Internal error, unspecified determinant")
      end if
      det = determinants%ground
      return
    end if

    select case(iDet)
    case(1)
      if (this%isGroundGuess) then
        det = determinants%ground
        return
      else if (this%isSpinPurify) then
        det = determinants%triplet
        return
      else
        det = determinants%mixed
        return
      end if
    case(2)
      if (this%isGroundGuess) then
        if (this%isSpinPurify) then
          det = determinants%triplet
          return
        else
          det = determinants%mixed
          return
        end if
      else
        if (this%isSpinPurify) then
          det = determinants%mixed
          return
        else
          call error("Internal error, unspecified determinant combination in Delta DFTB 2nd det")
        end if
      end if
    case(3)
      if (this%isGroundGuess .and. this%isSpinPurify) then
        det = determinants%mixed
        return
      else
        call error("Internal error, unspecified determinant combination in Delta DFTB 3rd det")
      end if
    case default
      call error("Internal error, unspecified determinant combination in Delta DFTB unknown det")
    end select

  end function whichDeterminant

  
  !> Spin Purifies Non-Aufbau excited state energy and forces
  subroutine postProcessDets(this, energy, derivs, tripletderivs, mixedderivs)

    !> Instance
    class(TDftbDeterminants), intent(inout) :: this

    !> energy components
    type(TEnergies), intent(inout) :: energy

    !> derivative components
    real(dp), intent(inout), optional:: derivs(:,:)

    !> Triplet state derivatives
    real(dp), intent(inout), optional :: tripletderivs(:,:)

    !> Spin contaminated derivatives
    real(dp), intent(inout), optional :: mixedderivs(:,:)

    this%isFinished = .true.

    if (.not.this%isNonAufbau) then
      return
    end if

    if (.not.this%isSpinPurify) then
      return
    end if

    ! Ziegler sum rule: E_S1 = 2E_mix - E_triplet
    energy%Etotal = 2.0_dp * energy%Emixed - energy%Etriplet
    energy%EMermin = 2.0_dp * energy%EmixMermin - energy%EtripMermin
    energy%Ezero = 2.0_dp * energy%EmixZero - energy%EtripZero
    energy%EGibbs = 2.0_dp * energy%EmixGibbs - energy%EtripGibbs
    energy%EForceRelated = 2.0_dp * energy%EmixForceRelated - energy%EtripForceRelated
    ! dE_S1 = 2dE_mix - dE_triplet
    if (present(derivs)) then
      derivs(:,:) = 2.0_dp * mixedDerivs - tripletDerivs
    endif

  end subroutine postProcessDets


  !> Initialised Time-independent excited state DFTB (TI-DFTB) conditions for determinant
  subroutine TDftbDeterminants_init(this, isNonAufbau, isSpinPurify, isGroundGuess, nEl)

    !> Instance
    type(TDftbDeterminants), intent(out) :: this

    !> Is this a non-aufbau filled TI calculation
    logical, intent(in) :: isNonAufbau

    !> Is this a spin purified calculation?
    logical, intent(in) :: isSpinPurify

    !> Should there be a ground state intial guess before Non-Aufbau calc?
    logical, intent(in) :: isGroundGuess

    !> Number of electrons in each spin channel
    real(dp), intent(in) :: nEl(:)

    this%isNonAufbau = isNonAufbau
    this%isGroundGuess = isGroundGuess
    this%isSpinPurify = isSpinPurify
    this%nEl = nEl

    this%isFinished = .false.

  end subroutine TDftbDeterminants_init


  !> Fillings for determinants
  subroutine detFilling(this, fillings, EBand, Ef, TS, E0, nElec, eigVals, tempElec, kWeights,&
      & iDistribFn)

    !> Instance
    class(TDftbDeterminants), intent(in) :: this

    !> Fillings (orbital, kpoint, spin)
    real(dp), intent(out) :: fillings(:,:,:)

    !> Band energies
    real(dp), intent(out) :: Eband(:)

    !> Fermi levels found for the given number of electrons on exit
    real(dp), intent(out) :: Ef(:)

    !> Band entropies
    real(dp), intent(out) :: TS(:)

    !> Band energies extrapolated to zero Kelvin
    real(dp), intent(out) :: E0(:)

    !> Number of electrons
    real(dp), intent(in) :: nElec(:)

    !> Eigenvalue of each level, kpoint and spin channel
    real(dp), intent(in) :: eigVals(:,:,:)

    !> Electronic temperature
    real(dp), intent(in) :: tempElec

    !> Weight of the k-points.
    real(dp), intent(in) :: kWeights(:)

    !> Selector for the distribution function
    integer, intent(in) :: iDistribFn

    real(dp), allocatable :: fillingsTmp(:,:,:,:), nElecFill(:)
    integer :: nSpinHams, nKPoints, nLevels, iS, iK, iConfig

    nLevels = size(fillings, dim=1)
    nKPoints = size(fillings, dim=2)
    nSpinHams = size(fillings, dim=3)

    nElecFill = nElec

    if (this%iDeterminant == determinants%mixed) then

      allocate(fillingsTmp(nLevels,nKPoints,2,3))
      fillingsTmp(:,:,:,:) = 0.0_dp

      do iConfig = 1, 3
        nElecFill(:2) = nElec(:2)
        select case(iConfig)
        case(1)
          nElecFill(1) = nElecFill(1) + 1.0_dp
        case(3)
          nElecFill(1) = nElecFill(1) - 1.0_dp
        end select
        ! Every spin channel (but not the k-points) filled up individually
        do iS = 1, nSpinHams
          call Efilling(Eband(iS:iS), Ef(iS), TS(iS:iS), E0(iS:iS),&
              & fillingsTmp(:,:,iS:iS, iConfig), eigvals(:,:,iS:iS), nElecFill(iS), tempElec,&
              & kWeights, iDistribFn)
        end do
      end do
      fillings(:,:,:) = fillingsTmp(:,:,:,1) - fillingsTmp(:,:,:,2) + fillingsTmp(:,:,:,3)

    else

      if (this%iDeterminant == determinants%triplet) then
        ! transfer an electron between spin channels
        nElecFill(1) = nElecFill(1) + 1.0_dp
        nElecFill(2) = nElecFill(2) - 1.0_dp
      end if

      ! Every spin channel (but not the k-points) filled up individually
      do iS = 1, nSpinHams
        call Efilling(Eband(iS:iS), Ef(iS), TS(iS:iS), E0(iS:iS), fillings(:,:,iS:iS),&
            & eigvals(:,:,iS:iS), nElecFill(iS), tempElec, kWeights, iDistribFn)
      end do

    end if

  end subroutine detFilling

end module dftbp_dftbdeterminants
