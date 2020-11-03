!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> main module for the DFTB+ API
module dftbp_mainapi
  use dftbp_accuracy, only : dp, mc
  use dftbp_assert
  use dftbp_coherence, only : checkExactCoherence, checkToleranceCoherence
  use dftbp_densedescr, only : TDenseDescr
  use dftbp_environment, only : TEnvironment
  use dftbp_initprogram
! use dftbp_initprogram, only : initProgramVariables, destructProgramVariables, TRefExtPot, 
!#:if WITH_SCALAPACK
!  use dftbp_initprogram, only : getDenseDescBlacs
!#:endif
  use dftbp_main, only : processGeometry
  use dftbp_message, only : error
  use dftbp_orbitals, only : TOrbitals
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy
#:if WITH_SCALAPACK
  use dftbp_scalapackfx, only : scalafx_getlocalshape
#:endif
  use dftbp_wrappedintr
  use dftbp_charmanip, only : newline
  implicit none
  private

! public :: initProgramVariables, destructProgramVariables
  public :: setGeometry, setQDepExtPotProxy, setExternalPotential, setExternalCharges
  public :: getEnergy, getGradients, getExtChargeGradients, getGrossCharges, getStressTensor
  public :: nrOfAtoms
  public :: updateDataDependentOnSpeciesOrdering, checkSpeciesNames

contains

  !> Sets up the atomic geometry
  subroutine setGeometry(env, globalData, coords, latVecs, coordOrigin)

    !> Instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TGlobalData), intent(inout) :: globalData

    !> atom coordinates
    real(dp), intent(in) :: coords(:,:)

    !> lattice vectors, if periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    !> coordinate origin, if periodic. If absent in that case, set to 0,0,0
    real(dp), intent(in), optional :: coordOrigin(:)

    ! Check data is consistent across MPI processes
  #:block DEBUG_CODE

    character(*), parameter :: routine = 'setGeometry'

    call checkToleranceCoherence(env, coords, "coords in "//routine, tol=1.e-10_dp)
    if (present(latVecs)) then
      call checkToleranceCoherence(env, latVecs, "latVecs in "//routine, tol=1.e-10_dp)
      call checkExactCoherence(env, globalData%tFracCoord, "tFracCoord in "//routine)
    endif
    if (present(coordOrigin)) then
      call checkToleranceCoherence(env, coordOrigin, "coordOrigin in "//routine, tol=1.e-10_dp)
    endif

  #:endblock DEBUG_CODE

    @:ASSERT(size(coords,1) == 3)

    globalData%coord0(:,:) = coords
    globalData%tCoordsChanged = .true.
    if (present(latVecs)) then
      globalData%latVec(:,:) = latVecs
      globalData%tLatticeChanged = .true.
      if (present(coordOrigin)) then
        globalData%origin = coordOrigin
      else
        globalData%origin = [0.0_dp,0.0_dp,0.0_dp]
      end if
    else
      globalData%tLatticeChanged = .false.
      @:ASSERT(.not.present(coordOrigin))
    end if

  end subroutine setGeometry


  !> Returns the free energy of the system at finite temperature
  subroutine getEnergy(env, globalData, merminEnergy)

    !> Instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TGlobalData), intent(in) :: globalData

    !> Resulting energy
    real(dp), intent(out) :: merminEnergy

    integer :: iDet

    call recalcGeometry(env, globalData)
    iDet = size(globalData%dftbEnergy)
    merminEnergy = globalData%dftbEnergy(iDet)%EMermin

  end subroutine getEnergy


  !> get forces on atoms
  subroutine getGradients(env, globalData, gradients)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TGlobalData), intent(inout) :: globalData

    !> resulting gradients wrt atom positions
    real(dp), intent(out) :: gradients(:,:)

    if (.not. globalData%tForces) then
      call error("Forces not available, you must initialise your calculator&
          & with forces enabled.")
    end if

    @:ASSERT(size(gradients,1) == 3)

    call recalcGeometry(env, globalData)
    gradients(:,:) = globalData%derivs

  end subroutine getGradients


  !> get stress tensor for unit cell
  subroutine getStressTensor(env, globalData, stress)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TGlobalData), intent(inout) :: globalData

    !> resulting gradients wrt atom positions
    real(dp), intent(out) :: stress(:,:)

    if (.not. globalData%tStress) then
      call error("Stress tensor not available, you must initialise your calculator with&
          & this property enabled.")
    end if

    call recalcGeometry(env, globalData)
    stress(:,:) = globalData%totalStress

  end subroutine getStressTensor


  !> get the gross (Mulliken projected) charges for atoms wrt neutral atoms
  subroutine getGrossCharges(env, globalData, atomCharges)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TGlobalData), intent(inout) :: globalData

    !> resulting charges
    real(dp), intent(out) :: atomCharges(:)

    call recalcGeometry(env, globalData)
    atomCharges(:) = sum(globalData%q0(:, :, 1) - globalData%qOutput(:, :, 1), dim=1)

    !> Pass to the charges of the excited state if relevant
    if (globalData%isLinResp) then
      atomCharges(:) = atomCharges(:) + globalData%dQAtomEx(:)
    end if

  end subroutine getGrossCharges


  !> Sets up an external population independent electrostatic potential.
  !>
  !> Sign convention: charge of electron is considered to be positive.
  !>
  subroutine setExternalPotential(globalData, atomPot, shellPot, potGrad)

    !> Instance
    type(TGlobalData), intent(inout) :: globalData

    !> Atomic external potential
    real(dp), intent(in), optional :: atomPot(:)

    !> Shell resolved electrostatic potential
    real(dp), intent(in), optional :: shellPot(:,:)

    !> Gradient of the electrostatic potential
    real(dp), intent(in), optional :: potGrad(:,:)

    ! Using explicit allocation instead of F2003 automatic ones in order to stop eventual
    ! shape mismatches already at this point rather than later deep in the main code
    if (present(atomPot)) then
      if (.not. allocated(globalData%refExtPot%atomPot)) then
        allocate(globalData%refExtPot%atomPot(globalData%nAtom, globalData%nSpin))
      end if
      @:ASSERT(all(shape(atomPot) == [globalData%nAtom]))
      globalData%refExtPot%atomPot(:,1) = atomPot
    end if
    if (present(shellPot)) then
      if (.not. allocated(globalData%refExtPot%shellPot)) then
        allocate(globalData%refExtPot%shellPot(globalData%orb%mShell, globalData%nAtom,&
            & globalData%nSpin))
      end if
      @:ASSERT(all(shape(shellPot) == [orb%mShell, globalData%nAtom]))
      globalData%refExtPot%shellPot(:,:,1) = shellPot
    end if
    if (present(potGrad)) then
      if (.not. allocated(globalData%refExtPot%potGrad)) then
        allocate(globalData%refExtPot%potGrad(3, globalData%nAtom))
      end if
      @:ASSERT(all(shape(potGrad) == [3, globalData%nAtom]))
      globalData%refExtPot%potGrad(:,:) = potGrad
    end if
    globalData%tExtField = .true.

  end subroutine setExternalPotential


  !> Sets up a generator for external population dependant potentials
  subroutine setQDepExtPotProxy(globalData, extPotProxy)

    !> Instance
    type(TGlobalData), intent(inout) :: globalData

    !> Generator for the external population dependant potential
    type(TQDepExtPotProxy), intent(in) :: extPotProxy

    globalData%qDepExtPot = extPotProxy

  end subroutine setQDepExtPotProxy


  !> Sets up external point charges
  subroutine setExternalCharges(globalData, chargeCoords, chargeQs, blurWidths)

    !> Instance
    type(TGlobalData), intent(inout) :: globalData

    !> Coordinates of the external charges
    real(dp), intent(in) :: chargeCoords(:,:)

    !> Charges of the external point charges (sign convention: electron is negative)
    real(dp), intent(in) :: chargeQs(:)

    !> Widths of the Gaussian for each charge used for blurring (0.0 = no blurring)
    real(dp), intent(in), optional :: blurWidths(:)

    globalData%tExtChrg = .true.
    if (globalData%tForces) then
      if (.not. allocated(globalData%chrgForces)) then
        allocate(globalData%chrgForces(3, size(chargeQs)))
      end if
    end if
    call globalData%sccCalc%setExternalCharges(chargeCoords, chargeQs, blurWidths=blurWidths)

  end subroutine setExternalCharges


  !> Returns the gradient acting on the external point charges
  subroutine getExtChargeGradients(globalData, chargeGradients)

    !> Instance
    type(TGlobalData), intent(in) :: globalData

    !> Gradients
    real(dp), intent(out) :: chargeGradients(:,:)

    @:ASSERT(globalData%tForces .and. allocated(globalData%chrgForces))

    chargeGradients(:,:) = globalData%chrgForces

  end subroutine getExtChargeGradients


  !> Obtains number of atoms in the system
  function nrOfAtoms(globalData)

    !> Instance
    type(TGlobalData), intent(in) :: globalData

    integer :: nrOfAtoms

    nrOfAtoms = globalData%nAtom

  end function nrOfAtoms


  !> Check that the order of speciesName remains constant Keeping speciesNames constant avoids the
  !> need to reset all of atomEigVal, referenceN0, speciesMass and SK parameters
  !>
  !> Even if nAtom is not conserved, it should be possible to know the total number of species
  !> types, nTypes, in a simulation and hence always keep speciesName constant
  function checkSpeciesNames(env, globalData, inputSpeciesName) result(tSpeciesNameChanged)

    !> dftb+ environment
    type(TEnvironment), intent(in) :: env

    !> Instance
    type(TGlobalData), intent(inout) :: globalData

    !> Labels of atomic species from external program
    character(len=*), intent(in) :: inputSpeciesName(:)

    !> Has speciesName changed?
    logical :: tSpeciesNameChanged

  #:block DEBUG_CODE

    ! Index over species
    integer :: ii

    ! Species index as string
    character(mc) :: i_str

    do ii = 1, size(inputSpeciesName)
      write(i_str, '(I0)') ii
      call checkExactCoherence(env, inputSpeciesName(ii), &
          & "inputSpeciesName element "//trim(i_str)//", in checkSpeciesNames")
    enddo

  #:endblock DEBUG_CODE

    tSpeciesNameChanged = any(globalData%speciesName /= inputSpeciesName)

  end function checkSpeciesNames


  !> When order of atoms changes, update arrays containing atom type indices,
  !> and all subsequent dependencies.
  !  Updated data returned via module use statements
  subroutine updateDataDependentOnSpeciesOrdering(env, globalData, inputSpecies)

    !> dftb+ environment
    type(TEnvironment), intent(in) :: env
 
    !> Instance
    type(TGlobalData), intent(inout) :: globalData

    !> types of the atoms (nAllAtom)
    integer, intent(in) :: inputSpecies(:)

    !> Dummy arguments. Won't be used if not allocated
    real(dp), allocatable :: initialCharges(:), initialSpins(:,:)
    type(TWrappedInt1), allocatable :: customOccAtoms(:)
    real(dp), allocatable :: customOccFillings(:,:)

    ! Check data is consistent across MPI processes
  #:block DEBUG_CODE

    character(*), parameter :: routine = 'updateDataDependentOnSpeciesOrdering'

    call checkExactCoherence(env, globalData%nAtom, "nAtom in "//routine)
    call checkExactCoherence(env, inputSpecies, "inputSpecies in "//routine)
    call checkExactCoherence(env, globalData%tSccCalc, "tSccCalc in" //routine)
    call checkExactCoherence(env, globalData%nSpin, "nSpin in "//routine)
    call checkExactCoherence(env, globalData%hamiltonianType, "hamiltonianType in "//routine)

  #:endblock DEBUG_CODE

    if(size(inputSpecies) /= globalData%nAtom)then
      call error("Number of atoms must be kept constant in simulation." // newline //&
          & "Instead call destruct and then fully re-initialize DFTB+.")
    endif

    globalData%species0 = inputSpecies
    globalData%mass =  updateAtomicMasses(globalData%species0)
    globalData%orb%nOrbAtom = updateAtomicOrbitals(globalData%species0)

    ! if atom species change, dense matrix indexing needs updating
    call globalData%getDenseDescCommon(globalData%orb, globalData%nAtom, globalData%t2Component,&
        & globalData%denseDesc)

    ! Used in partial charge initialisation
    call globalData%setEquivalencyRelations(globalData%species0, globalData%sccCalc,&
        & globalData%orb, globalData%onSiteElements, globalData%iEqOrbitals,&
        & globalData%iEqBlockDFTBU, globalData%iEqBlockOnSite, globalData%iEqBlockDFTBULS,&
        & globalData%iEqBlockOnSiteLS, globalData%nIneqOrb, globalData%nMixElements)
  #:if WITH_SCALAPACK
    call updateBLACSDecomposition(env, globalData, globalData%denseDesc)
    call reallocateHSArrays(env, globalData%denseDesc, globalData%HSqrCplx, globalData%SSqrCplx,&
        & globalData%eigVecsCplx, globalData%HSqrReal, globalData%SSqrReal, globalData%eigVecsReal)
  #:endif

    ! If atomic order changes, partial charges need to be initialised, else wrong charge will be
    ! associated with each atom
    call globalData%initializeReferenceCharges(globalData%species0, globalData%referenceN0,&
        & globalData%orb, customOccAtoms, customOccFillings, globalData%q0, globalData%qshell0)
    call globalData%setNElectrons(globalData%q0, globalData%nrChrg, globalData%nrSpinPol,&
        & globalData%nEl, globalData%nEl0)
    call globalData%initializeCharges(globalData%species0, globalData%speciesName, globalData%orb,&
        & globalData%nEl, globalData%iEqOrbitals, globalData%nIneqOrb, globalData%nMixElements,&
        & initialSpins, initialCharges, globalData%nrChrg, globalData%q0, globalData%qInput,&
        & globalData%qOutput, globalData%qDiff, globalData%qInpRed, globalData%qOutRed,&
        & globalData%qDiffRed, globalData%qBlockIn, globalData%qBlockOut, globalData%qiBlockIn,&
        & globalData%qiBlockOut)

  end subroutine updateDataDependentOnSpeciesOrdering


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update order of nr. atomic orbitals for each atom, orb%nOrbAtom
  function updateAtomicOrbitals(globalData, species0) result(nOrbAtomReordered)

    !> Instance
    type(TGlobalData), intent(in) :: globalData

    !> Type of the atoms (nAtom)
    integer, intent(in)  :: species0(:)

    !> Nr. of orbitals for each atom (nAtom)
    integer, allocatable :: nOrbAtomReordered(:)

    @:ASSERT(globalData%nAtom == size(species0))
    allocate(nOrbAtomReordered(globalData%nAtom))
    nOrbAtomReordered(:) = globalData%orb%nOrbSpecies(species0(:))

  end function updateAtomicOrbitals


  !> Update atomic masses
  function updateAtomicMasses(globalData, species0) result(massReordered)

    !> Instance
    type(TGlobalData), intent(in) :: globalData

    !> Type of the atoms (nAtom)
    integer, intent(in)  :: species0(:)

    !> List of atomic masses (nAtom)
    real(dp), allocatable :: massReordered(:)

    @:ASSERT(size(globalData%speciesMass) == maxval(species0))
    allocate(massReordered(globalData%nAtom))
    massReordered = globalData%speciesMass(species0)

  end function updateAtomicMasses

#:if WITH_SCALAPACK

  !> Update dense matrix descriptor for H and S in BLACS decomposition
  subroutine updateBLACSDecomposition(env, globalData, denseDesc)

    !> Environment settings
    type(TEnvironment), intent(in)    :: env

    !> Instance
    type(TGlobalData), intent(inout) :: globalData

    !> Dense matrix descriptor for H and S
    type(TDenseDescr), intent(inout) :: denseDesc


    ! Specificaly, denseDesc uses orb%nOrbAtom
    call globalData%getDenseDescCommon(globalData%orb, globalData%nAtom, globalData%t2Component,&
        & denseDesc)
    call globalData%getDenseDescBlacs(env, env%blacs%rowBlockSize, env%blacs%columnBlockSize,&
        & denseDesc)

  end subroutine updateBLACSDecomposition


  !> Reassign Hamiltonian, overlap and eigenvector arrays
  !
  ! If nAtom is constant and one is running without BLACS, this should not be required
  ! hence preprocessed out
  ! May require extending if ((nAtom not constant) and (not BLACS))
  subroutine reallocateHSArrays(env, globalData, denseDesc, HSqrCplx, SSqrCplx, eigVecsCplx,&
      & HSqrReal, SSqrReal ,eigVecsReal)

    !> Environment instance
    type(TEnvironment), intent(in) :: env

    !> Instance
    type(TGlobalData), intent(in) :: globalData

    !> Dense matrix descriptor for H and S
    type(TDenseDescr), intent(in) :: denseDesc

    !> Square dense hamiltonian storage for cases with k-points
    complex(dp), allocatable, intent(inout) :: HSqrCplx(:,:)

    !> Square dense overlap storage for cases with k-points
    complex(dp), allocatable, intent(inout) :: SSqrCplx(:,:)

    !> Complex eigenvectors
    complex(dp), allocatable, intent(inout) :: eigvecsCplx(:,:,:)

    !> Square dense hamiltonian storage
    real(dp), allocatable, intent(inout) :: HSqrReal(:,:)

    !> Square dense overlap storage
    real(dp), allocatable, intent(inout) :: SSqrReal(:,:)

    !> Real eigenvectors
    real(dp), allocatable, intent(inout) :: eigvecsReal(:,:,:)

    ! Local variables
    integer :: nLocalRows, nLocalCols, nLocalKS

    !Retrieved from index array for spin and k-point index
    nLocalKS = size(globalData%parallelKS%localKS, dim=2)

    !Get nLocalRows and nLocalCols
    call scalafx_getlocalshape(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, nLocalRows, nLocalCols)

    if (globalData%t2Component .or. .not. globalData%tRealHS) then
      ! Complex
      if( (size(HSqrCplx,1) /= nLocalRows) .or. (size(HSqrCplx,2) /= nLocalCols) )then
        deallocate(HSqrCplx)
        deallocate(SSqrCplx)
        deallocate(eigVecsCplx)
        allocate(HSqrCplx(nLocalRows, nLocalCols))
        allocate(SSqrCplx(nLocalRows, nLocalCols))
        allocate(eigVecsCplx(nLocalRows, nLocalCols, nLocalKS))
      endif
    else
      ! Real
      if( (size(HSqrReal,1) /= nLocalRows) .or. (size(HSqrReal,2) /= nLocalCols) )then
        deallocate(HSqrReal)
        deallocate(SSqrReal)
        deallocate(eigVecsReal)
        allocate(HSqrReal(nLocalRows, nLocalCols))
        allocate(SSqrReal(nLocalRows, nLocalCols))
        allocate(eigVecsReal(nLocalRows, nLocalCols, nLocalKS))
      endif

    end if

  end subroutine reallocateHSArrays

#:endif

  !> re-evaluate the energy/forces if the geometry changes
  subroutine recalcGeometry(env, globalData)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TGlobalData), intent(inout) :: globalData

    logical :: tStopScc, tExitGeoOpt

    if (globalData%tLatticeChanged .or. globalData%tCoordsChanged) then
      call processGeometry(env, globalData, 1, 1, .false., tStopScc, tExitGeoOpt)
      globalData%tLatticeChanged = .false.
      globalData%tCoordsChanged = .false.
    end if

  end subroutine recalcGeometry

end module dftbp_mainapi
