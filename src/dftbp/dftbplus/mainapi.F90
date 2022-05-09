!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> main module for the DFTB+ API
module dftbp_dftbplus_mainapi
  use dftbp_common_accuracy, only : dp, mc
  use dftbp_common_coherence, only : checkExactCoherence, checkToleranceCoherence
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_status, only : TStatus
  use dftbp_dftbplus_initprogram, only : TDftbPlusMain, initReferenceCharges, initElectronNumbers
  use dftbp_dftbplus_main, only : processGeometry
  use dftbp_dftbplus_qdepextpotproxy, only : TQDepExtPotProxy
  use dftbp_io_charmanip, only : newline
  use dftbp_io_message, only : error
  use dftbp_timedep_timeprop, only : initializeDynamics, doTdStep
  use dftbp_type_densedescr, only : TDenseDescr
  use dftbp_type_orbitals, only : TOrbitals
  use dftbp_type_wrappedintr, only : TWrappedInt1
#:if WITH_SCALAPACK
  use dftbp_dftbplus_initprogram, only : getDenseDescBlacs
  use dftbp_extlibs_scalapackfx, only : scalafx_getlocalshape
#:endif
  implicit none

  private
  public :: setGeometry, setQDepExtPotProxy, setExternalPotential, setExternalCharges
  public :: getEnergy, getGradients, getExtChargeGradients, getGrossCharges, getCM5Charges
  public :: getElStatPotential, getStressTensor, nrOfAtoms, nrOfKPoints, getAtomicMasses
  public :: updateDataDependentOnSpeciesOrdering, checkSpeciesNames
  public :: initializeTimeProp, doOneTdStep, setTdElectricField, setTdCoordsAndVelos, getTdForces


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


  !> get the CM5 charges
  subroutine getCM5Charges(env, main, atomCharges)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> resulting charges
    real(dp), intent(out) :: atomCharges(:)

    !> number of neighbours for all atoms
    integer, allocatable :: nNeigh(:)

    !> handle the case that CM5 was not added in the input
    if (.not. allocated(main%cm5Cont)) then
      call error("CM5 analysis has not been carried out.")
    end if

    call recalcGeometry(env, main)
    if (.not. allocated(main%cm5Cont%cm5)) then
      call error("CM5 could not be calculated.")
    end if
    atomCharges(:) = sum(main%q0(:, :, 1) - main%qOutput(:, :, 1), dim=1) + main%cm5Cont%cm5

    !> Pass to the charges of the excited state if relevant
    if (main%isLinResp) then
      atomCharges(:) = atomCharges + main%dQAtomEx
    end if

  end subroutine getCM5Charges


  !>  get electrostatic potential at specified points
  subroutine getElStatPotential(env, main, pot, locations)

    !> instance
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> Resulting potentials
    real(dp), intent(out) :: pot(:)

    !> Sites to calculate potential
    real(dp), intent(in) :: locations(:,:)

    !> Default potential softening
    real(dp) :: epsSoften = 1E-6

    call main%scc%getInternalElStatPotential(pot, env, locations, epsSoften)

  end subroutine getElStatPotential


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
    main%isExtField = .true.

    ! work around for lack (at the moment) for a flag to re-calculate ground state even if
    ! geometries are unchanged.
    main%tCoordsChanged = .true.

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
    call main%scc%setExternalCharges(chargeCoords, chargeQs, blurWidths=blurWidths)

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


  !> Obtains number of k-points in the system (1 if not a repeating structure)
  function nrOfKPoints(main)

    !> Instance
    type(TDftbPlusMain), intent(in) :: main

    integer :: nrOfKPoints

    nrOfKPoints = main%nKPoint

  end function nrOfKPoints


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

    if (main%atomOrderMatters) then
      call error("This DftbPlus instance can not cope with atom reordering (by initialization)")
    end if

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
    call initReferenceCharges(main%species0, main%orb, main%referenceN0, main%nSpin, main%q0,&
        & main%qShell0)
    call initElectronNumbers(main%q0, main%nrChrg, main%nrSpinPol, main%nSpin, main%orb,&
        & main%nEl0, main%nEl)
    call main%initializeCharges()

  end subroutine updateDataDependentOnSpeciesOrdering


  !> After calculation of the ground state, this subroutine initializes the variables
  !> and the initial step of the propagators for electron and nuclear dynamics
  subroutine initializeTimeProp(env, main, dt, tdFieldThroughAPI, tdCoordsAndVelosThroughAPI)

    !> dftb+ environment
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> time step
    real(dp), intent(in) :: dt

    !> field will be provided through the API?
    logical, intent(in) :: tdFieldThroughAPI

    !> coords and velocities will be provided at each step through the API?
    logical, intent(in) :: tdCoordsAndVelosThroughAPI

    type(TStatus) :: errStatus

    if (allocated(main%electronDynamics)) then
      main%electronDynamics%tdFieldThroughAPI = tdFieldThroughAPI
      if (tdCoordsAndVelosThroughAPI) then
        if (main%electronDynamics%tIons) then
          main%electronDynamics%tdCoordsAndVelosThroughAPI = tdCoordsAndVelosThroughAPI
        else
          call error("Setting coordinates and velocities at each step is allowed only for&
              & simulations with ion dynamics enabled")
        end if
      end if

      main%electronDynamics%dt = dt
      main%electronDynamics%iCall = 1
      call initializeDynamics(main%electronDynamics, main%boundaryCond, main%coord0, main%orb,&
          & main%neighbourList, main%nNeighbourSK, main%denseDesc%iAtomStart, main%iSparseStart,&
          & main%img2CentCell, main%skHamCont, main%skOverCont, main%ints, env, main%coord,&
          & main%H0, main%spinW, main%tDualSpinOrbit, main%xi, main%thirdOrd, main%dftbU,&
          & main%onSiteElements, main%refExtPot, main%solvation, main%eFieldScaling, main%rangeSep,&
          & main%referenceN0, main%q0, main%repulsive, main%iAtInCentralRegion, main%eigvecsReal,&
          & main%eigvecsCplx, main%filling, main%qDepExtPot, main%tFixEf, main%Ef, main%latVec,&
          & main%invLatVec, main%iCellVec, main%rCellVec, main%cellVec, main%species,&
          & main%electronicSolver, errStatus)
      if (errStatus%hasError()) then
        call error(errStatus%message)
      end if
    else
      call error("Electron dynamics not enabled, please initialize the calculator&
          & including the ElectronDynamics block")
    end if

  end subroutine initializeTimeProp


  !> After calling initializeTimeProp, this subroutine performs one timestep of
  !> electron and nuclear (if IonDynamics enabled) dynamics.
  subroutine doOneTdStep(env, main, iStep, dipole, energy, atomNetCharges,&
      & coordOut, force, occ)

    !> dftb+ environment
    type(TEnvironment), intent(inout) :: env

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> present step of dynamics
    integer, intent(in) :: iStep

    !> Dipole moment
    real(dp), optional, intent(out) :: dipole(:,:)

    !> total energy
    real(dp), optional, intent(out) :: energy

    !> Negative gross charge
    real(dp), optional, intent(out) :: atomNetCharges(:,:)

    !> atomic coordinates
    real(dp), optional, intent(out) :: coordOut(:,:)

    !> forces (3, nAtom)
    real(dp), optional, intent(out) :: force(:,:)

    !> molecular orbital projected populations
    real(dp), optional, intent(out) :: occ(:)

    type(TStatus) :: errStatus

    if (main%electronDynamics%tPropagatorsInitialized) then
      call doTdStep(main%electronDynamics, main%boundaryCond, iStep, main%coord0, main%orb,&
          & main%neighbourList, main%nNeighbourSK,main%denseDesc%iAtomStart, main%iSparseStart,&
          & main%img2CentCell, main%skHamCont, main%skOverCont, main%ints, env, main%coord,&
          & main%q0, main%referenceN0, main%spinW, main%tDualSpinOrbit, main%xi, main%thirdOrd,&
          & main%dftbU, main%onSiteElements, main%refExtPot, main%solvation, main%eFieldScaling,&
          & main%rangeSep, main%repulsive, main%iAtInCentralRegion, main%tFixEf, main%Ef,&
          & main%electronicSolver, main%qDepExtPot, errStatus)

      if (errStatus%hasError()) then
        call error(errStatus%message)
      end if

      if (present(dipole)) then
        dipole(:,:) = main%electronDynamics%dipole
      end if
      if (present(energy)) then
        energy = main%electronDynamics%energy%eSCC
      end if
      if (present(atomNetCharges)) then
        atomNetCharges(:,:) = main%electronDynamics%deltaQ
      end if
      if (present(coordOut)) then
        coordOut(:,:) = main%coord0
      end if
      if (present(force)) then
        force(:,:) = main%electronDynamics%totalForce
      end if
      if (present(occ)) then
        occ(:) = main%electronDynamics%occ
      end if
    else
      call error("Propagators for dynamics not initialize, please call initializeTimeProp()&
          & first.")
    end if

  end subroutine doOneTdStep


  !> sets electric field for td propagation
  subroutine setTdElectricField(main, field)

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    ! electric field components
    real(dp), intent(in) :: field(3)

    main%electronDynamics%presentField(:) = field
    main%electronDynamics%tdFieldIsSet = .true.

  end subroutine setTdElectricField


  !> sets coordinates and velos for td propagation
  subroutine setTdCoordsAndVelos(main, coords, velos)

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    ! coordinates
    real(dp), intent(in) :: coords(3, main%nAtom)

    ! velocities
    real(dp), intent(in) :: velos(3, main%nAtom)

    main%electronDynamics%coordNew(:,:) = coords
    main%electronDynamics%movedVelo(:,:) = velos(:, main%electronDynamics%indMovedAtom)
    main%electronDynamics%tdCoordsAndVelosAreSet = .true.

  end subroutine setTdCoordsAndVelos


  !> gets atomic forces from time dependent propagation
  subroutine getTdForces(main, forces)

    !> Instance
    type(TDftbPlusMain), intent(inout) :: main

    !> forces (3, nAtom)
    real(dp), intent(out) :: forces(:,:)

    forces(:,:) = main%electronDynamics%totalForce

  end subroutine getTdForces


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
        & main%denseDesc, main%isSparseReorderRequired)

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

    !> Status of operation
    type(TStatus) :: errStatus

    if (main%tLatticeChanged .or. main%tCoordsChanged) then
      call processGeometry(main, env, 1, 1, .false., tStopScc, tExitGeoOpt, errStatus)
      if (errStatus%hasError()) then
        call error(errStatus%message)
      end if
      main%tLatticeChanged = .false.
      main%tCoordsChanged = .false.
    end if

  end subroutine recalcGeometry

end module dftbp_dftbplus_mainapi
