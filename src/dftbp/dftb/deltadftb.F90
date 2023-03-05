!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#! Template parameters for array dimensioning
#:set TEMPLATE_PARAMS = [('0', ''), ('1', '(:)'), ('3', '(:,:,:)'), ('4', '(:,:,:,:)')]

!> Routines for (time independent excited) TI-DFTB
module dftbp_dftb_deltadets
  use dftbp_common_accuracy, only : dp
  use dftbp_dftb_energytypes, only : TEnergies
  use dftbp_dftb_etemp, only : Efilling
  use dftbp_io_message, only : error
  implicit none

  private
  public :: deltaDeterminants, TDeltaDeterminants, TDeltaDeterminants_init


  !> Name space for Delta-SCF determinants
  type :: TDeltaDeterminantsEnum

    !> Conventional ground state
    integer :: ground = 0

    !> Triplet configuration
    integer :: triplet = 1

    !> S1 state, spin contaminated
    integer :: mixedS1 = 2

    !> S2 state, electron in higher level, spin contaminated
    integer :: mixedS2e = 3

    !> S2 state, hole in deeper level, spin contaminated
    integer :: mixedS2h = 4

  end type TDeltaDeterminantsEnum


  !> Actual values for Delta-SCF determinants
  type(TDeltaDeterminantsEnum), parameter :: deltaDeterminants = TDeltaDeterminantsEnum()


  !> Control type for Delta-SCF / TI-excitations
  type TDeltaDeterminants

    !> Is this a non-Aufbau filled state
    logical :: isNonAufbau

    !> Current determinant being solved
    integer :: iDet = 0

    !> Should the non-Aufbau fillings be spin purified
    logical :: isSpinPurify

    !> Should a ground state case be calculated as well
    logical :: isGroundGuess

    !> Has the calculation finished and results are now ready to use
    logical :: isFinished

    !> which determinant holds the ground state (0 if none)
    integer :: iGround

    !> which determinant holds the triplet state (0 if none)
    integer :: iTriplet

    !> which determinant holds the (spin contaminated) excited state
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

  end type TDeltaDeterminants


  !> Ziegler rule for purification
  interface zieglerSum
    #:for SUFFIX, _ in TEMPLATE_PARAMS
    module procedure applyZiegler${SUFFIX}$
    #:endfor
  end interface zieglerSum


contains

  !> Counts number of determinants to be evaluated for a property calculation
  pure function nDeterminant(this) result(nDet)

    !> Instance
    class(TDeltaDeterminants), intent(in) :: this

    integer :: nDet

    nDet = size(this%determinants)

  end function nDeterminant


  !> Converts determinant number into what type it should be
  function whichDeterminant(this) result(det)

    !> Instance
    class(TDeltaDeterminants), intent(in) :: this

    integer :: det

    det = this%determinants(this%iDet)

  end function whichDeterminant


  !> Spin purifies Non-Aufbau excited state properties
  subroutine postProcessDets(this, energies, qOutput, qDets, qBlockOut, qBlockDets, dipoleMoment,&
      & stress, tripletStress, mixedStress, derivs, tripletderivs, mixedderivs)

    !> Instance
    class(TDeltaDeterminants), intent(inout) :: this

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
          call zieglerSum(eM%Etotal, eT%Etotal, eM1%Etotal)
          call zieglerSum(eM%Ezero, eT%Ezero, eM1%Ezero)
          call zieglerSum(eM%EMermin, eT%EMermin, eM1%EMermin)
          call zieglerSum(eM%EGibbs, eT%EGibbs, eM1%EGibbs)
          call zieglerSum(eM%EForceRelated, eT%EForceRelated, eM1%EForceRelated)
          call zieglerSum(eM%TS, eT%TS, eM1%TS)
          call zieglerSum(eM%EBand, eT%EBand, eM1%EBand)
          call zieglerSum(eM%E0, eT%E0, eM1%E0)
          call zieglerSum(eM%atomRep, eT%atomRep, eM1%atomRep)
          call zieglerSum(eM%atomNonScc, eT%atomNonScc, eM1%atomNonScc)
          call zieglerSum(eM%atomScc, eT%atomScc, eM1%atomScc)
          call zieglerSum(eM%atomSpin, eT%atomSpin, eM1%atomSpin)
          call zieglerSum(eM%atomLS, eT%atomLS, eM1%atomLS)
          call zieglerSum(eM%atomDftbU, eT%atomDftbU, eM1%atomDftbU)
          call zieglerSum(eM%atomExt, eT%atomExt, eM1%atomExt)
          call zieglerSum(eM%atomElec, eT%atomElec, eM1%atomElec)
          call zieglerSum(eM%atomDisp, eT%atomDisp, eM1%atomDisp)
          call zieglerSum(eM%atomOnSite, eT%atomOnSite, eM1%atomOnSite)
          call zieglerSum(eM%atomHalogenX, eT%atomHalogenX, eM1%atomHalogenX)
          call zieglerSum(eM%atom3rd, eT%atom3rd, eM1%atom3rd)
          call zieglerSum(eM%atomSolv, eT%atomSolv, eM1%atomSolv)
          call zieglerSum(eM%atomTotal, eT%atomTotal, eM1%atomTotal)
        end associate
      end if

      if (allocated(qDets)) then
        if (this%isSpinPurify) then
          ! S1 = 2 mix - triplet
          qOutput(:,:,:) = 2.0_dp * qDets(:,:,:,this%iMixed) - qDets(:,:,:,this%iTriplet)
        else
          ! this is copying from the last calculated determinant at the moment, so is redundant
          qOutput(:,:,:) = qDets(:,:,:,this%iMixed)
        end if
      end if

      if (allocated(qBlockDets)) then
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
        else
          dipoleMoment(:,this%iFinal) = dipoleMoment(:,this%iMixed)
        end if
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

    end if

  end subroutine postProcessDets


  !> Initialised Time-independent excited state conditions for determinant
  subroutine TDeltaDeterminants_init(this, isNonAufbau, isSpinPurify, isGroundGuess)

    !> Instance
    type(TDeltaDeterminants), intent(out) :: this

    !> Is this a non-aufbau filled TI calculation
    logical, intent(in) :: isNonAufbau

    !> Is this a spin purified calculation?
    logical, intent(in) :: isSpinPurify

    !> Should there be a ground state initial guess before Non-Aufbau calc?
    logical, intent(in) :: isGroundGuess

    if (.not.isNonAufbau .and. (isGroundGuess .or. isSpinPurify) ) then
      call error("Delta DFTB internal error - setting request without non-Aufbau fillings")
    end if

    this%isNonAufbau = isNonAufbau
    this%isGroundGuess = isGroundGuess
    this%isSpinPurify = isSpinPurify

    ! set to zero as initially unused.
    this%iGround = 0
    this%iTriplet = 0
    this%iMixed = 0

    ! assume first determinant
    this%iDet = 1

    this%determinants = [integer ::]

    if (isNonAufbau) then
      if (isGroundGuess) then
        ! first determinant
        this%iGround = 1
        this%determinants = [deltaDeterminants%ground]
      end if
      if (isSpinPurify) then
        ! if required, then its after the ground state (if that is requested)
        this%iTriplet = size(this%determinants) + 1
        this%determinants = [this%determinants, deltaDeterminants%triplet]
      end if
      ! last determinant
      this%iMixed = size(this%determinants) + 1
      this%determinants = [this%determinants, deltaDeterminants%mixedS1]

      if (isSpinPurify) then
        ! need one extra element to store the resulting purified results
        this%iFinal = size(this%determinants) + 1
      else
        ! only storing the different determinant results
        this%iFinal = size(this%determinants)
      end if
    else
      this%iGround = 1
      this%determinants = [deltaDeterminants%ground]
      this%iFinal = size(this%determinants)
    end if

    this%isFinished = .false.

  end subroutine TDeltaDeterminants_init


  !> Fillings for determinants in Delta SCF
  subroutine detFilling(this, fillings, EBand, Ef, TS, E0, nElec, eigVals, tempElec, kWeights,&
      & iDistribFn)

    !> Instance
    class(TDeltaDeterminants), intent(in) :: this

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

    if (this%whichDeterminant() == deltaDeterminants%mixedS1) then

      allocate(fillingsTmp(nLevels, nKPoints, 2, 3))
      fillingsTmp(:,:,:,:) = 0.0_dp

      if (.true.) then
        do iConfig = 1, 3
          nElecFill(1) = nElec(1)
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
      else
      end if
      fillings(:,:,:) = fillingsTmp(:,:,:,1) - fillingsTmp(:,:,:,2) + fillingsTmp(:,:,:,3)

    else

      if (this%whichDeterminant() == deltaDeterminants%triplet) then
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


  #:for SUFFIX, DIMS in TEMPLATE_PARAMS

  !> Apply Zieglers rule for ${SUFFIX}$-dim arrays.
  pure subroutine applyZiegler${SUFFIX}$(eM, eT, eF)

    !> Mixed determinant value
    real(dp), intent(in) :: eM${DIMS}$

    !> Triplet determinant value
    real(dp), intent(in) :: eT${DIMS}$

    ! Final resulting value
    real(dp), intent(out) :: eF${DIMS}$

    eF${DIMS}$ = 2.0_dp * eM - eT

  end subroutine applyZiegler${SUFFIX}$

  #:endfor

end module dftbp_dftb_deltadets
