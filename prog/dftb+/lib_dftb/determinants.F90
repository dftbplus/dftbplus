!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Routines for (time independent excited) TI-DFTB
module dftbp_determinants
  use dftbp_accuracy, only : dp
  use dftbp_energytypes, only : TEnergies
  use dftbp_message, only : error
  use dftbp_etemp, only : Efilling
  use dftbp_assert
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

    !> which determinant holds the ground state (0 if none)
    integer :: iGround

    !> which determinant holds the triplet state (0 if none)
    integer :: iTriplet

    !> which determinant holds the (spin contaminated) S1 state
    integer :: iMixed

    !> Resulting final determinant
    integer :: iFinal

    !> List of determinants to be calculated
    integer, allocatable :: determinants(:)

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

    nDet = size(this%determinants)

  end function nDeterminant


  !> Converts determinant number into what type it should be
  function whichDeterminant(this, iDet) result(det)

    !> Instance
    class(TDftbDeterminants), intent(in) :: this

    !> Number of current determinant
    integer, intent(in) :: iDet

    integer :: det

    if (iDet > size(this%determinants)) then
      call error("Internal error: invalid determinant")
    endif

    det = this%determinants(iDet)

  end function whichDeterminant


  !> Spin Purifies Non-Aufbau excited state energy and forces
  subroutine postProcessDets(this, energies, qOutput, qDets, qBlockOut, qBlockDets, dipoleMoment,&
      & stress, tripletStress, mixedStress, derivs, tripletderivs, mixedderivs)

    !> Instance
    class(TDftbDeterminants), intent(inout) :: this

    !> energy components for whatever determinants are present
    type(TEnergies), intent(inout) :: energies(:)

    !> Charges
    real(dp), intent(inout), allocatable :: qOutput(:,:,:)

    !> Block charges from determinants
    real(dp), intent(in), allocatable :: qDets(:,:,:,:)

    !> Block charges
    real(dp), intent(inout), allocatable :: qBlockOut(:,:,:,:)

    !> Charges from determinants
    real(dp), intent(in), allocatable :: qBlockDets(:,:,:,:,:)

    !> dipole moment
    real(dp), intent(inout), allocatable :: dipoleMoment(:,:)

    !> stress tensor
    real(dp), intent(inout) :: stress(:,:)

    !> Triplet stress
    real(dp), intent(inout), optional :: tripletStress(:,:)

    !> Spin contaminated stress
    real(dp), intent(inout), optional :: mixedStress(:,:)

    !> derivatives for atom positions
    real(dp), intent(inout), optional:: derivs(:,:)

    !> Triplet state derivatives
    real(dp), intent(inout), optional :: tripletDerivs(:,:)

    !> Spin contaminated derivatives
    real(dp), intent(inout), optional :: mixedDerivs(:,:)

    this%isFinished = .true.

    if (this%isNonAufbau) then
      if (this%isSpinPurify) then
        associate(eM => energies(this%iMixed), eM1 => energies(this%iMixed + 1),&
            & eT => energies(this%iTriplet))
          call applyZiegler(eM%Etotal, eT%Etotal, eM1%Etotal)
          call applyZiegler(eM%Ezero, eT%Ezero, eM1%Ezero)
          call applyZiegler(eM%EMermin, eT%EMermin, eM1%EMermin)
          call applyZiegler(eM%EGibbs, eT%EGibbs, eM1%EGibbs)
          call applyZiegler(eM%EForceRelated, eT%EForceRelated, eM1%EForceRelated)
          call applyZieglerAlloc(eM%TS, eT%TS, eM1%TS)
          call applyZieglerAlloc(eM%EBand, eT%EBand, eM1%EBand)
          call applyZieglerAlloc(eM%E0, eT%E0, eM1%E0)
          call applyZieglerAlloc(eM%atomRep, eT%atomRep, eM1%atomRep)
          call applyZieglerAlloc(eM%atomNonScc, eT%atomNonScc, eM1%atomNonScc)
          call applyZieglerAlloc(eM%atomScc, eT%atomScc, eM1%atomScc)
          call applyZieglerAlloc(eM%atomSpin, eT%atomSpin, eM1%atomSpin)
          call applyZieglerAlloc(eM%atomLS, eT%atomLS, eM1%atomLS)
          call applyZieglerAlloc(eM%atomDftbU, eT%atomDftbU, eM1%atomDftbU)
          call applyZieglerAlloc(eM%atomExt, eT%atomExt, eM1%atomExt)
          call applyZieglerAlloc(eM%atomElec, eT%atomElec, eM1%atomElec)
          call applyZieglerAlloc(eM%atomDisp, eT%atomDisp, eM1%atomDisp)
          call applyZieglerAlloc(eM%atomOnSite, eT%atomOnSite, eM1%atomOnSite)
          call applyZieglerAlloc(eM%atomHalogenX, eT%atomHalogenX, eM1%atomHalogenX)
          call applyZieglerAlloc(eM%atom3rd, eT%atom3rd, eM1%atom3rd)
          call applyZieglerAlloc(eM%atomSolv, eT%atomSolv, eM1%atomSolv)
          call applyZieglerAlloc(eM%atomTotal, eT%atomTotal, eM1%atomTotal)
        end associate
      end if
    end if

    if (allocated(qDets)) then
      @:ASSERT(this%isNonAufbau)
      if (this%isSpinPurify) then
        ! S1 = 2 mix - triplet
        qOutput(:,:,:) = 2.0_dp * qDets(:,:,:,this%iMixed) - qDets(:,:,:,this%iTriplet)
      else
        ! this is copying from the last calculated determinant at the moment, so is redundant
        qOutput(:,:,:) = qDets(:,:,:,this%iMixed)
      end if
    end if

    if (allocated(qBlockDets)) then
      @:ASSERT(this%isNonAufbau)
      if (this%isSpinPurify) then
        ! S1 = 2 mix - triplet
        qBlockOut(:,:,:,:) = 2.0_dp*qBlockDets(:,:,:,:,this%iMixed)&
            & - qBlockDets(:,:,:,:,this%iTriplet)
      else
        ! this is copying from the last calculated determinant at the moment, so is redundant
        qBlockOut(:,:,:,:) = qBlockDets(:,:,:,:,this%iMixed)
      end if
    end if

    if (allocated(dipoleMoment)) then
      if (this%isSpinPurify) then
        dipoleMoment(:,this%iFinal) = 2.0_dp * dipoleMoment(:,this%iMixed)&
            & - dipoleMoment(:,this%iTriplet)
      end if
    end if

    if (.not.this%isNonAufbau) then
      return
    end if

    if (present(derivs)) then
      if (this%isSpinPurify) then
        ! dE_S1 = 2dE_mix - dE_triplet
        derivs(:,:) = 2.0_dp * mixedDerivs - tripletDerivs
      else
        derivs(:,:) = mixedDerivs
      end if
    endif

    if (present(mixedStress)) then
      if (this%isSpinPurify) then
        ! dE_S1 = 2dE_mix - dE_triplet
        stress(:,:) = 2.0_dp * mixedStress - tripletStress
      else
        derivs(:,:) = mixedStress
      end if
    end if

  end subroutine postProcessDets


  !> Initialised Time-independent excited state DFTB (TI-DFTB) conditions for determinant
  subroutine TDftbDeterminants_init(this, isNonAufbau, isSpinPurify, isGroundGuess, nEl, dftbEnergy)

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

    !> Energy terms for each determinant (and total if post processing something other than ground
    !> state)
    type(TEnergies), allocatable, intent(out) :: dftbEnergy(:)

    if (.not.isNonAufbau .and. (isGroundGuess .or. isSpinPurify) ) then
      call error("Delta DFTB internal error - setting request without non-Aufbau fillings")
    end if

    this%isNonAufbau = isNonAufbau
    this%isGroundGuess = isGroundGuess
    this%isSpinPurify = isSpinPurify

    this%nEl = nEl

    ! set to zero as initially unused.
    this%iGround = 0
    this%iTriplet = 0
    this%iMixed = 0

    ! assume first determinant
    this%iDeterminant = 1

    this%determinants = [integer ::]

    if (isNonAufbau) then
      if (isGroundGuess) then
        ! first determinant
        this%iGround = 1
        this%determinants = [determinants%ground]
      end if
      if (isSpinPurify) then
        ! if required, then its after the ground state (if that is requested)
        this%iTriplet = size(this%determinants) + 1
        this%determinants = [this%determinants, determinants%triplet]
      end if
      ! last determinant
      this%iMixed = size(this%determinants) + 1
      this%determinants = [this%determinants, determinants%mixed]

      if (isSpinPurify) then
        ! need one extra element to store the resulting purified results
        allocate(dftbEnergy(size(this%determinants) + 1))
      else
        ! only storing the different determinant results
        allocate(dftbEnergy(size(this%determinants)))
      end if
    else
      this%iGround = 1
      this%determinants = [determinants%ground]
      allocate(dftbEnergy(size(this%determinants)))
    end if

    this%iFinal = size(dftbEnergy)

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
    integer :: nSpinHams, nKPoints, nLevels, iS, iConfig

    nLevels = size(fillings, dim=1)
    nKPoints = size(fillings, dim=2)
    nSpinHams = size(fillings, dim=3)

    nElecFill = nElec

    if (this%whichDeterminant(this%iDeterminant) == determinants%mixed) then
      allocate(fillingsTmp(nLevels, nKPoints, 2, 3))
      fillingsTmp(:,:,:,:) = 0.0_dp

      do iConfig = 1, 3
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
        nElecFill(1) = nElec(1)
      end do
      fillings(:,:,:) = fillingsTmp(:,:,:,1) - fillingsTmp(:,:,:,2) + fillingsTmp(:,:,:,3)

    else
      if (this%whichDeterminant(this%iDeterminant) == determinants%triplet) then
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


  !> Apply Zieglers rule.
  elemental subroutine applyZiegler(eM, eT, eF)

    !> Mixed energy value
    real(dp), intent(in) :: eM

    !> Triple energy value
    real(dp), intent(in) :: eT

    ! Final resulting energy value
    real(dp), intent(out) :: eF

    eF = 2.0_dp * eM - eT

  end subroutine applyZiegler


  !> Apply Zieglers rule for 1D arrays if they are allocated
  pure subroutine applyZieglerAlloc(eM, eT, eF)

    !> Mixed energy value
    real(dp), intent(in), optional :: eM(:)

    !> Triple energy value
    real(dp), intent(in), optional :: eT(:)

    !> Final resulting energy value
    real(dp), intent(out), optional :: eF(:)

    if (present(eF)) then
      call applyZiegler(eM, eT, eF)
    end if

  end subroutine applyZieglerAlloc


end module dftbp_determinants
