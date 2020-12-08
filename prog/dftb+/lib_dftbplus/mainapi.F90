!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> main module for the DFTB+ API
module dftbp_mainapi
  use dftbp_assert
  use dftbp_accuracy, only : dp, mc
  use dftbp_coherence, only : checkExactCoherence, checkToleranceCoherence
  use dftbp_densedescr, only : TDenseDescr
  use dftbp_environment, only : TEnvironment
  use dftbp_initprogram, only : TDftbPlusMain
#:if WITH_SCALAPACK
  use dftbp_initprogram, only : getDenseDescBlacs
#:endif
  use dftbp_main, only : processGeometry
  use dftbp_message, only : error
  use dftbp_orbitals, only : TOrbitals
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy
#:if WITH_SCALAPACK
  use dftbp_scalapackfx, only : scalafx_getlocalshape
#:endif
  use dftbp_wrappedintr, only : TWrappedInt1
  use dftbp_charmanip, only : newline
  implicit none
  private

  public :: setGeometry, setQDepExtPotProxy, setExternalPotential, setExternalCharges
  public :: getEnergy, getGradients, getExtChargeGradients, getGrossCharges, getStressTensor
  public :: nrOfAtoms, getAtomicMasses
  public :: updateDataDependentOnSpeciesOrdering, checkSpeciesNames


contains

  !> Sets up the atomic geometry
  subroutine setGeometry(env, main, coords, latVecs, coordOrigin)

    !> Instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

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
      call checkExactCoherence(env, main%tFracCoord, "tFracCoord in "//routine)
    endif
    if (present(coordOrigin)) then
      call checkToleranceCoherence(env, coordOrigin, "coordOrigin in "//routine, tol=1.e-10_dp)
    endif

  #:endblock DEBUG_CODE

    @:ASSERT(size(coords,1) == 3)

    main%coord0(:,:) = coords
    main%tCoordsChanged = .true.
    if (present(latVecs)) then
      main%latVec(:,:) = latVecs
      main%tLatticeChanged = .true.
      if (present(coordOrigin)) then
        main%origin = coordOrigin
      else
        main%origin = [0.0_dp,0.0_dp,0.0_dp]
      end if
    else
      main%tLatticeChanged = .false.
      @:ASSERT(.not.present(coordOrigin))
    end if

  end subroutine setGeometry


  !> Returns the free energy of the system at finite temperature
  subroutine getEnergy(env, main, merminEnergy)

    !> Instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> Resulting energy
    real(dp), intent(out) :: merminEnergy

    integer :: iDet

    call recalcGeometry(env, main)
    iDet = size(main%dftbEnergy)
    merminEnergy = main%dftbEnergy(iDet)%EMermin

  end subroutine getEnergy


  !> get forces on atoms
  subroutine getGradients(env, main, gradients)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> resulting gradients wrt atom positions
    real(dp), intent(out) :: gradients(:,:)

    if (.not. main%tForces) then
      call error("Forces not available, you must initialise your calculator&
          & with forces enabled.")
    end if

    @:ASSERT(size(gradients,1) == 3)

    call recalcGeometry(env, main)
    gradients(:,:) = main%derivs

  end subroutine getGradients


  !> get stress tensor for unit cell
  subroutine getStressTensor(env, main, stress)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> resulting gradients wrt atom positions
    real(dp), intent(out) :: stress(:,:)

    if (.not. main%tStress) then
      call error("Stress tensor not available, you must initialise your calculator with&
          & this property enabled.")
    end if

    call recalcGeometry(env, main)
    stress(:,:) = main%totalStress

  end subroutine getStressTensor


  !> get the gross (Mulliken projected) charges for atoms wrt neutral atoms
  subroutine getGrossCharges(env, main, atomCharges)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> resulting charges
    real(dp), intent(out) :: atomCharges(:)

    call recalcGeometry(env, main)
    atomCharges(:) = sum(main%q0(:, :, 1) - main%qOutput(:, :, 1), dim=1)

    !> Pass to the charges of the excited state if relevant
    if (main%isLinResp) then
      atomCharges(:) = atomCharges(:) + main%dQAtomEx(:)
    end if

  end subroutine getGrossCharges


  !> Sets up an external population independent electrostatic potential.
  !>
  !> Sign convention: charge of electron is considered to be positive.
  !>
  subroutine setExternalPotential(main, atomPot, shellPot, potGrad)

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> Atomic external potential
    real(dp), intent(in), optional :: atomPot(:)

    !> Shell resolved electrostatic potential
    real(dp), intent(in), optional :: shellPot(:,:)

    !> Gradient of the electrostatic potential
    real(dp), intent(in), optional :: potGrad(:,:)

    ! Using explicit allocation instead of F2003 automatic ones in order to stop eventual
    ! shape mismatches already at this point rather than later deep in the main code
    if (present(atomPot)) then
      if (.not. allocated(main%refExtPot%atomPot)) then
        allocate(main%refExtPot%atomPot(main%nAtom, main%nSpin))
      end if
      @:ASSERT(all(shape(atomPot) == [main%nAtom]))
      main%refExtPot%atomPot(:,1) = atomPot
    end if
    if (present(shellPot)) then
      if (.not. allocated(main%refExtPot%shellPot)) then
        allocate(main%refExtPot%shellPot(main%orb%mShell, main%nAtom,&
            & main%nSpin))
      end if
      @:ASSERT(all(shape(shellPot) == [main%orb%mShell, main%nAtom]))
      main%refExtPot%shellPot(:,:,1) = shellPot
    end if
    if (present(potGrad)) then
      if (.not. allocated(main%refExtPot%potGrad)) then
        allocate(main%refExtPot%potGrad(3, main%nAtom))
      end if
      @:ASSERT(all(shape(potGrad) == [3, main%nAtom]))
      main%refExtPot%potGrad(:,:) = potGrad
    end if
    main%tExtField = .true.

  end subroutine setExternalPotential


  !> Sets up a generator for external population dependant potentials
  subroutine setQDepExtPotProxy(main, extPotProxy)

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> Generator for the external population dependant potential
    type(TQDepExtPotProxy), intent(in) :: extPotProxy

    main%qDepExtPot = extPotProxy

  end subroutine setQDepExtPotProxy


  !> Sets up external point charges
  subroutine setExternalCharges(main, chargeCoords, chargeQs, blurWidths)

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> Coordinates of the external charges
    real(dp), intent(in) :: chargeCoords(:,:)

    !> Charges of the external point charges (sign convention: electron is negative)
    real(dp), intent(in) :: chargeQs(:)

    !> Widths of the Gaussian for each charge used for blurring (0.0 = no blurring)
    real(dp), intent(in), optional :: blurWidths(:)

    main%tExtChrg = .true.
    if (main%tForces) then
      if (.not. allocated(main%chrgForces)) then
        allocate(main%chrgForces(3, size(chargeQs)))
      end if
    end if
    call main%sccCalc%setExternalCharges(chargeCoords, chargeQs, blurWidths=blurWidths)

  end subroutine setExternalCharges


  !> Returns the gradient acting on the external point charges
  subroutine getExtChargeGradients(main, chargeGradients)

    !> Instance
    type(TDftbPlusMain), intent(in) :: main

    !> Gradients
    real(dp), intent(out) :: chargeGradients(:,:)

    @:ASSERT(main%tForces .and. allocated(main%chrgForces))

    chargeGradients(:,:) = main%chrgForces

  end subroutine getExtChargeGradients


  !> Obtains number of atoms in the system
  function nrOfAtoms(main)

    !> Instance
    type(TDftbPlusMain), intent(in) :: main

    integer :: nrOfAtoms

    nrOfAtoms = main%nAtom

  end function nrOfAtoms


  !> Check that the order of speciesName remains constant Keeping speciesNames constant avoids the
  !> need to reset all of atomEigVal, referenceN0, speciesMass and SK parameters
  !>
  !> Even if nAtom is not conserved, it should be possible to know the total number of species
  !> types, nTypes, in a simulation and hence always keep speciesName constant
  function checkSpeciesNames(env, main, inputSpeciesName) result(tSpeciesNameChanged)

    !> dftb+ environment
    type(TEnvironment), intent(in) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

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

    tSpeciesNameChanged = any(main%speciesName /= inputSpeciesName)

  end function checkSpeciesNames


  !> When order of atoms changes, update arrays containing atom type indices,
  !> and all subsequent dependencies.
  !  Updated data returned via module use statements
  subroutine updateDataDependentOnSpeciesOrdering(env, main, inputSpecies)

    !> dftb+ environment
    type(TEnvironment), intent(in) :: env
 
    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> types of the atoms (nAllAtom)
    integer, intent(in) :: inputSpecies(:)

    !> Dummy arguments. Won't be used if not allocated
    real(dp), allocatable :: initialCharges(:), initialSpins(:,:)
    type(TWrappedInt1), allocatable :: customOccAtoms(:)
    real(dp), allocatable :: customOccFillings(:,:)

    ! Check data is consistent across MPI processes
  #:block DEBUG_CODE

    character(*), parameter :: routine = 'updateDataDependentOnSpeciesOrdering'

    call checkExactCoherence(env, main%nAtom, "nAtom in "//routine)
    call checkExactCoherence(env, inputSpecies, "inputSpecies in "//routine)
    call checkExactCoherence(env, main%tSccCalc, "tSccCalc in" //routine)
    call checkExactCoherence(env, main%nSpin, "nSpin in "//routine)
    call checkExactCoherence(env, main%hamiltonianType, "hamiltonianType in "//routine)

  #:endblock DEBUG_CODE

    if(size(inputSpecies) /= main%nAtom)then
      call error("Number of atoms must be kept constant in simulation." // newline //&
          & "Instead call destruct and then fully re-initialize DFTB+.")
    endif

    main%species0 = inputSpecies
    main%mass = updateAtomicMasses(main)
    main%orb%nOrbAtom = updateAtomicOrbitals(main)

    ! if atom species change, dense matrix indexing needs updating
    call main%getDenseDescCommon()

    ! Used in partial charge initialisation
    call main%setEquivalencyRelations()
  #:if WITH_SCALAPACK
    call updateBLACSDecomposition(env, main)
    call reallocateHSArrays(env, main, main%denseDesc, main%HSqrCplx,&
        & main%SSqrCplx, main%eigVecsCplx, main%HSqrReal, main%SSqrReal,&
        & main%eigVecsReal)
  #:endif

    ! If atomic order changes, partial charges need to be initialised, else wrong charge will be
    ! associated with each atom
    call main%initializeReferenceCharges(customOccAtoms, customOccFillings)
    call main%setNElectrons()
    call main%initializeCharges(initialSpins, initialCharges)

  end subroutine updateDataDependentOnSpeciesOrdering


  !> Obtains mass for each atom in the system
  subroutine getAtomicMasses(main, outMass)

    !> Instance
    type(TDftbPlusMain), intent(in) :: main

    real(dp), intent(out) :: outMass(main%nAtom)

    outMass = main%mass

  end subroutine getAtomicMasses



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update order of nr. atomic orbitals for each atom, orb%nOrbAtom
  function updateAtomicOrbitals(main) result(nOrbAtomReordered)

    !> Instance
    type(TDftbPlusMain), intent(in) :: main

    !> Nr. of orbitals for each atom (nAtom)
    integer, allocatable :: nOrbAtomReordered(:)

    @:ASSERT(main%nAtom == size(main%species0))
    allocate(nOrbAtomReordered(main%nAtom))
    nOrbAtomReordered(:) = main%orb%nOrbSpecies(main%species0(:))

  end function updateAtomicOrbitals


  !> Update atomic masses
  function updateAtomicMasses(main) result(massReordered)

    !> Instance
    type(TDftbPlusMain), intent(in) :: main

    !> List of atomic masses (nAtom)
    real(dp), allocatable :: massReordered(:)

    @:ASSERT(size(main%speciesMass) == maxval(main%species0))
    allocate(massReordered(main%nAtom))
    massReordered = main%speciesMass(main%species0)

  end function updateAtomicMasses

#:if WITH_SCALAPACK

  !> Update dense matrix descriptor for H and S in BLACS decomposition
  subroutine updateBLACSDecomposition(env, main)

    !> Environment settings
    type(TEnvironment), intent(in)    :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main


    ! Specificaly, denseDesc uses orb%nOrbAtom
    call main%getDenseDescCommon()
    call getDenseDescBlacs(env, env%blacs%rowBlockSize, env%blacs%columnBlockSize,&
        & main%denseDesc)

  end subroutine updateBLACSDecomposition


  !> Reassign Hamiltonian, overlap and eigenvector arrays
  !
  ! If nAtom is constant and one is running without BLACS, this should not be required
  ! hence preprocessed out
  ! May require extending if ((nAtom not constant) and (not BLACS))
  subroutine reallocateHSArrays(env, main, denseDesc, HSqrCplx, SSqrCplx, eigVecsCplx,&
      & HSqrReal, SSqrReal, eigVecsReal)

    !> Environment instance
    type(TEnvironment), intent(in) :: env

    !> Instance
    type(TDftbPlusMain), intent(in) :: main

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
    nLocalKS = size(main%parallelKS%localKS, dim=2)

    !Get nLocalRows and nLocalCols
    call scalafx_getlocalshape(env%blacs%orbitalGrid, denseDesc%blacsOrbSqr, nLocalRows, nLocalCols)

    if (main%t2Component .or. .not. main%tRealHS) then
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
  subroutine recalcGeometry(env, main)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    logical :: tStopScc, tExitGeoOpt

    if (main%tLatticeChanged .or. main%tCoordsChanged) then
      call processGeometry(main, env, 1, 1, .false., tStopScc, tExitGeoOpt)
      main%tLatticeChanged = .false.
      main%tCoordsChanged = .false.
    end if

  end subroutine recalcGeometry

end module dftbp_mainapi
