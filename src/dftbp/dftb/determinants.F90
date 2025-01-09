!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Routines for (time independent excited) TI-DFTB
module dftbp_dftb_determinants
  use dftbp_common_accuracy, only : dp
  use dftbp_common_globalenv, only : stdOut, withMpi
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_energytypes, only : TEnergies
  use dftbp_dftb_etemp, only : Efilling
  use dftbp_io_message, only : error
  use dftbp_dftb_periodic, only: TNeighbourList
  use dftbp_type_commontypes, only: TOrbitals
  use dftbp_type_densedescr, only: TDenseDescr
  use dftbp_common_environment, only : TEnvironment
  use dftbp_math_lapackroutines, only : gesvd
  use dftbp_math_matrixops, only : adjointLowerTriangle
  use dftbp_dftb_populations, only : mulliken
  use dftbp_dftb_sparse2dense, only : packHS
  implicit none

  private
  public :: determinants, TDftbDeterminants_init, TDftbDeterminants
  public :: tiTDM, tiTraDip


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


  !> Names of the determinants, matching TDeterminantsEnum
  character(len=2), parameter :: detNames(0:2) = ['s0', 't1', 's1']


  !> Control type for Delta DFTB / TI-DFTB
  type TDftbDeterminants

    !> Is this a non-Aufbau filling
    logical :: isNonAufbau = .false.

    !> Should the non-Aufbau fillings be spin purified
    logical :: isSpinPurify

    !> Should a ground state case be calculated as well
    logical :: isGroundGuess

    !> Should the transition dipole moment between the ground and excited state be computed?
    logical :: isTDM

    !> Current determinant being solved
    integer :: iDeterminant

    !> Number of electrons in each spin channel
    real(dp), allocatable :: nEl(:)

    !> Has the calculation finished and results are now ready to use
    logical :: isFinished

    !> Which determinant holds the ground state (0 if none)
    integer :: iGround

    !> Which determinant holds the triplet state (0 if none)
    integer :: iTriplet

    !> Which determinant holds the (spin contaminated) S1 state
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
    procedure :: determinantName

  end type TDftbDeterminants


contains

  !> Counts number of determinants to be evaluated
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

    !> Determinant number corresponding to detNames / determinants variables in this class
    integer :: det

    if (iDet > size(this%determinants)) then
      call error("Internal error: invalid determinant")
    endif

    det = this%determinants(iDet)

  end function whichDeterminant


  !> Converts determinant number into what it is called
  function determinantName(this, iDet) result(det)

    !> Instance
    class(TDftbDeterminants), intent(in) :: this

    !> Number of current determinant
    integer, intent(in) :: iDet

    !> Determinant name corresponding to detNames in this class
    character(2) :: det

    if (iDet > size(this%determinants)) then
      call error("Internal error: invalid determinant")
    endif

    det = detNames(this%determinants(iDet))

  end function determinantName


  !> Spin Purifies Non-Aufbau excited state energy and forces
  subroutine postProcessDets(this, energies, qOutput, qDets, qBlockOut, qBlockDets, dipoleMoment,&
      & transitionDipoleMoment, groundFill, mixedFill, stress, tripletStress, mixedStress, derivs,&
      & tripletderivs, mixedderivs)

    !> Instance
    class(TDftbDeterminants), intent(inout) :: this

    !> Energy components for whatever determinants are present
    type(TEnergies), intent(inout) :: energies(:)

    !> Charges
    real(dp), intent(inout), allocatable :: qOutput(:,:,:)

    !> Block charges from determinants
    real(dp), intent(in), allocatable :: qDets(:,:,:,:)

    !> Block charges
    real(dp), intent(inout), allocatable :: qBlockOut(:,:,:,:)

    !> Charges from determinants
    real(dp), intent(in), allocatable :: qBlockDets(:,:,:,:,:)

    !> Dipole moment
    real(dp), intent(inout), allocatable :: dipoleMoment(:,:)

    !> Transition dipole moment
    real(dp), intent(inout), allocatable :: transitionDipoleMoment(:)

    !> Fillings for ground state (TI-DFTB transition dipoles)
    real(dp), intent(inout), allocatable :: groundFill(:,:,:)

    !> Fillings for mixed-spin state (TI-DFTB transition dipoles)
    real(dp), intent(inout), allocatable :: mixedFill(:,:,:)

    !> Stress tensor
    real(dp), intent(inout) :: stress(:,:)

    !> Triplet stress
    real(dp), intent(inout), optional :: tripletStress(:,:)

    !> Spin contaminated stress
    real(dp), intent(inout), optional :: mixedStress(:,:)

    !> Derivatives for atom positions
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
  subroutine TDftbDeterminants_init(this, isNonAufbau, isSpinPurify, isGroundGuess, isTDM, nEl,&
      & dftbEnergy, errStatus)

    !> Instance
    type(TDftbDeterminants), intent(out) :: this

    !> Is this a non-aufbau filled TI calculation
    logical, intent(in) :: isNonAufbau

    !> Is this a spin purified calculation?
    logical, intent(in) :: isSpinPurify

    !> Should there be a ground state initial guess before Non-Aufbau calc?
    logical, intent(in) :: isGroundGuess

    !> Should the transition dipole moment between the ground and excited state be computed?
    logical, intent(in) :: isTDM

    !> Number of electrons in each spin channel
    real(dp), intent(in) :: nEl(:)

    !> Energy terms for each determinant (and total if post processing something other than ground
    !> state)
    type(TEnergies), allocatable, intent(out) :: dftbEnergy(:)

    type(TStatus), intent(out) :: errStatus

    if (.not.isNonAufbau .and. (isGroundGuess .or. isSpinPurify) ) then
      @:RAISE_ERROR(errStatus, -1, "Delta DFTB internal error - setting request without non-Aufbau&
          & fillings")
    end if

    this%isNonAufbau = isNonAufbau
    this%isGroundGuess = isGroundGuess
    this%isSpinPurify = isSpinPurify
    this%isTDM = isTDM
    if (withMpi .and. this%isTDM) then
      @:RAISE_ERROR(errStatus, -1, "Delta DFTB transition dipole not yet available for MPI enabled&
          & builds")
    end if

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
        ! if required, then it is after the ground state (if that is requested)
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

      allocate(fillingsTmp(nLevels, nKPoints, 2, 3), source=0.0_dp)

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


  !> Apply Ziegler's spin purification rule.
  elemental subroutine applyZiegler(eM, eT, eF)

    !> Mixed energy value
    real(dp), intent(in) :: eM

    !> Triple energy value
    real(dp), intent(in) :: eT

    !> Final resulting energy value
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


  !> Evaluate transition density matrix for time independent excitation
  subroutine tiTDM(groundC, excitedC, transitionDM, groundFill, mixedFill)

    !> DFTB ground state MO eigenvectors
    real(dp), intent(inout) :: groundC(:,:)

    !> TI-DFTB mixed-determinant MO eigenvectors
    real(dp), intent(inout) :: excitedC(:,:)

    !> Transition density matrix in corresponding orbital basis
    real(dp), intent(inout) :: transitionDM(:,:)

    !> DFTB ground state occupation numbers
    real(dp), intent(in) :: groundFill(:)

    !> TI-DFTB mixed-determinant occupation numbers
    real(dp), intent(in) :: mixedFill(:)

    !! Left singular vectors (corresponding orbitals)
    real(dp), allocatable :: u(:,:)

    !! Right (transposed) singular vectors (corresponding orbitals)
    real(dp), allocatable :: vt(:,:)

    !! Temporary storage for matrix products in this subroutine
    real(dp), allocatable :: work(:,:)

    !! Singular values
    real(dp), allocatable :: sigma(:)

    !! DFTB ground state MO eigenvectors in transformed (CO) basis
    real(dp), allocatable :: groundMOs(:,:)

    !! TI-DFTB excited state MO eigenvectors in transformed (CO) basis
    real(dp), allocatable :: excitedMOs(:,:)

    integer :: jj, nElec, nOrb, nGrndMOs, nExMOs

    nOrb = size(excitedC,dim=1)
    @:ASSERT(nOrb == size(groundC,dim=1))
    nGrndMOs = size(groundC, dim=2)
    nExMOs = size(excitedC, dim=2)

    allocate(u(nOrb, nExMOs))
    allocate(vt(nOrb, nExMOs))
    allocate(work(nOrb, nExMOs))
    allocate(sigma(nOrb))
    allocate(excitedMOs(nOrb, nExMOs))
    allocate(groundMOs(nOrb, nGrndMOs))

    !! Determine how many orbitals with non-negligible occupations to track (assumes Fermi
    !! occupation)
    do jj = size(groundFill), 1, -1
      nElec = jj
      if (abs(groundFill(jj)) >= epsilon(1.0_dp)) then
        exit
      end if
    end do

    !! Compute singular value decomposition of non-orthogonal MO-MO overlap
    work = matmul(transpose(groundC),excitedC)
    call gesvd(work,u,sigma,vt)

    !! Rotate ground MOs into corresponding orbital basis and isolate the occupied MOs
    work = matmul(groundC, u)
    groundMOs(:,:) = 0.0_dp
    do jj = 1, nElec
      if (abs(groundFill(jj)) >= epsilon(1.0_dp)) then
        groundMOs(:,jj) = groundFill(jj) * work(:,jj)
      end if
    end do

    !! Rotate excited MOs into corresponding orbital basis and isolate the occupied MOs
    work = matmul(excitedC,transpose(vt))
    excitedMOs(:,:) = 0.0_dp
    do jj = 1, nElec
      excitedMOs(:,jj) = 0.0_dp
      if (abs(mixedFill(jj)) >= epsilon(1.0_dp)) then
        excitedMOs(:,jj) = mixedFill(jj) * work(:,jj)
      end if
    end do
    if (mixedFill(nElec) == 0.0_dp) then
      do jj = nElec, size(mixedFill) ! Could just as well start from nElec+1 here (see condition)
        if (abs(mixedFill(jj)) >= epsilon(1.0_dp)) then
          excitedMOs(:,nElec) = mixedFill(jj) * work(:,jj)
        end if
      end do
    end if

    transitionDM = matmul(groundMOs,transpose(excitedMOs))

  end subroutine tiTDM


  !> Evaluate transition dipole for time independent excitation
  subroutine tiTraDip(tiMatG, tiMatE, tiMatPT, neighbourlist, nNeighbourSK, orb, denseDesc,&
      & iSparseStart, img2CentCell, rhoPrimSize, over, tiTraCharges, transitionDipoleMoment, q0,&
      & coord, iAtInCentralRegion, groundFill, mixedFill, env)

    !> DFTB ground state MO eigenvectors
    real(dp), intent(inout) :: tiMatG(:,:,:)

    !> TI-DFTB mixed-determinant MO eigenvectors
    real(dp), intent(inout) :: tiMatE(:,:,:)

    !> Transition density matrix in corresponding orbital basis
    real(dp), intent(inout) :: tiMatPT(:,:)

    !> List of neighbouring atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours of each real atom
    integer, intent(in) :: nNeighbourSK(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Dense matrix descriptor
    type(TDenseDescr), intent(in) :: denseDesc

    !> Index array for the start of atomic blocks in sparse arrays
    integer, intent(in) :: iSparseStart(:,:)

    !> Image atoms to their equivalent in the central cell
    integer, allocatable, intent(in) :: img2CentCell(:)

    !> Sparse density matrix storage(used only for size)
    real(dp), allocatable, intent(in) :: rhoPrimSize(:,:)

    !> sparse overlap matrix
    real(dp), intent(in) :: over(:)

    !> Transition Charges
    real(dp), intent(inout) :: tiTraCharges(:,:)

    ! > Transition dipole moment
    real(dp), intent(inout) :: transitionDipoleMoment(:)

    !> reference atomic charges
    real(dp), intent(in) :: q0(:,:,:)

    !> atomic coordinates
    real(dp), intent(in) :: coord(:,:)

    !> atoms in the central cell (or device region if transport)
    integer, intent(in) :: iAtInCentralRegion(:)

    !> ground state filling for TDM
    real(dp), intent(in) :: groundFill(:,:,:)

    !> mixed-spin excited state filling for TDM
    real(dp), intent(in) :: mixedFill(:,:,:)

    !> Sparse representation for transition density matrix
    real(dp), allocatable :: rhoPrim(:,:)

    !> Dense representation for transition density matrix
    real(dp), allocatable :: tiTransitionDensity(:,:,:)

    !> Center of atomic transition charges
    real(dp), allocatable :: tiTrCenter(:)

    !> Sum of transition charges
    real(dp) :: sumOfTraCharges

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    integer :: nAtom, ii, iAtom, ia, na, jj, m, iSpin

    na = size(tiMatE, dim=1)
    nAtom = size(orb%nOrbAtom)

    allocate(rhoPrim(size(rhoPrimSize, dim=1), size(rhoPrimSize, dim=2)))
    allocate(tiTransitionDensity(size(tiMatE, dim=1), size(tiMatE, dim=2), size(tiMatE, dim=3)))
    allocate(tiTrCenter(size(coord, dim=1)))

    ! Reset transition dipole moment
    transitionDipoleMoment(:) = 0.0_dp

    ! Reset transition charges (Mulliken pops for transition density matrix)
    tiTraCharges(:,:) = 0.0_dp

    ! Corresponding orbital transformation of stored ground and excited MOs
    do iSpin = 1, size(groundFill, dim=3)

      call tiTDM(tiMatG(:,:,iSpin), tiMatE(:,:,iSpin), tiMatPT, groundFill(:,1,iSpin),&
          & mixedFill(:,1,iSpin))

      ! Matrix of transition density
      tiTransitionDensity(:,:, iSpin) = tiMatPT
      call adjointLowerTriangle(tiTransitionDensity(:,:, iSpin))

      ! Pack TI-DFTB transition density matrix (only lower triangle referenced)
      rhoPrim(:,iSpin) = 0.0_dp
      call packHS(rhoPrim(:,iSpin), tiTransitionDensity(:,:,iSpin), neighbourlist%iNeighbour,&
          & nNeighbourSK, orb%mOrb, denseDesc%iAtomStart, iSparseStart, img2CentCell)

      call mulliken(env, tiTraCharges, over, rhoPrim(:,iSpin), orb, neighbourlist%iNeighbour,&
          & nNeighbourSK, img2CentCell, iSparseStart)

    #:block DEBUG_CODE
      ! Write out the transition charges (for debugging)
      write(stdOut,*) "QM Transition Charges for spin ", iSpin
      do ia = 1, size(tiTraCharges, dim=2)
        write(stdOut, *) tiTraCharges(:, ia)
      end do
    #:endblock DEBUG_CODE

      ! Shift center of transition charges to address gauge dependence
      sumOfTraCharges = 0.0_dp
      tiTrCenter(:) = 0.0_dp
      do ii = 1, size(iAtInCentralRegion)
        iAtom = iAtInCentralRegion(ii)
        sumOfTraCharges = sumOfTraCharges + sum(tiTraCharges(:,iAtom))
      end do
      do ii = 1, size(iAtInCentralRegion)
        iAtom = iAtInCentralRegion(ii)
        tiTrCenter(:) = tiTrCenter(:) + sum(tiTraCharges(:,iAtom))*coord(:,iAtom)&
            &/sumOfTraCharges
      end do

      do ii = 1, size(iAtInCentralRegion)
        iAtom = iAtInCentralRegion(ii)
        transitionDipoleMoment(:) = transitionDipoleMoment(:)&
            & + sum(q0(:, iAtom, iSpin) - tiTraCharges(:, iAtom))&
            & * (coord(:, iAtom)-tiTrCenter(:))
      end do

    end do

  end subroutine tiTraDip

end module dftbp_dftb_determinants
