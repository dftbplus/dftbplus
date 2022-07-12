!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides DFTB+ API for MM-type high level access
module dftbp_mmapi
  use iso_fortran_env, only : output_unit
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment, TEnvironment_init
  use dftbp_common_globalenv, only : initGlobalEnv, destructGlobalEnv, instanceSafeBuild, withMpi
  use dftbp_dftbplus_hsdhelpers, only : doPostParseJobs
  use dftbp_dftbplus_initprogram, only: TDftbPlusMain
  use dftbp_dftbplus_inputdata, only : TInputData
  use dftbp_dftbplus_mainapi, only : doOneTdStep, checkSpeciesNames, nrOfAtoms, nrOfKPoints,&
      & setExternalPotential, getTdForces, setTdCoordsAndVelos, setTdElectricField,&
      & initializeTimeProp, updateDataDependentOnSpeciesOrdering, getAtomicMasses,&
      & getGrossCharges, getCM5Charges, getElStatPotential, getExtChargeGradients, getStressTensor,&
      & getGradients, getEnergy, setQDepExtPotProxy, setExternalCharges, setGeometry
  use dftbp_dftbplus_parser, only : TParserFlags, rootTag, parseHsdTree, readHsdFile
  use dftbp_dftbplus_qdepextpotgen, only : TQDepExtPotGen, TQDepExtPotGenWrapper
  use dftbp_dftbplus_qdepextpotproxy, only : TQDepExtPotProxy, TQDepExtPotProxy_init
  use dftbp_extlibs_xmlf90, only : fnode, createDocumentNode, createElement, appendChild
  use dftbp_io_charmanip, only : newline
  use dftbp_io_hsdutils, only : getChild
  use dftbp_io_message, only: error
  use dftbp_type_linkedlist, only : TListString, init, append, get, len, asArray
  use dftbp_type_typegeometry, only : TGeometry
  implicit none
  private

  public :: TDftbPlus, getDftbPlusBuild, getDftbPlusApi
  public :: TDftbPlus_init, TDftbPlus_destruct
  public :: TDftbPlusAtomList
  public :: TDftbPlusInput
  public :: TQDepExtPotGen
  public :: getMaxAngFromSlakoFile, convertAtomTypesToSpecies


  !> list of QM atoms and species for DFTB+ calculation
  type :: TDftbPlusAtomList
    !> number of atoms
    integer :: nAtom
    !> linked list of chemical symbols of elements (species names), size=nSpecies
    type(TListString) :: speciesNames
    !> array of species for each atom, size=nAtom
    integer, allocatable :: species(:)
  contains
    !> read list of atoms
    procedure :: get => TDftbPlusAtomList_get
    !> insert the list of atoms into the input data structure
    procedure :: add => TDftbPlusAtomList_addToInpData
  end type TDftbPlusAtomList


  !> input tree for DFTB+ calculation
  type :: TDftbPlusInput
    !> tree for HSD format input
    type(fnode), pointer :: hsdTree => null()
  contains
    !> obtain the root of the tree of input
    procedure :: getRootNode => TDftbPlusInput_getRootNode
  end type TDftbPlusInput


  !> A DFTB+ calculation
  type :: TDftbPlus
    private
    type(TEnvironment), allocatable :: env
    type(TDftbPlusMain), allocatable :: main
    logical :: tInit = .false.
  contains
    !> read input from a file
    procedure :: getInputFromFile => TDftbPlus_getInputFromFile
    !> get an empty input to populate
    procedure :: getEmptyInput => TDftbPlus_getEmptyInput
    !> set up a DFTB+ calculator from input tree
    procedure :: setupCalculator => TDftbPlus_setupCalculator
    !> set/replace the geometry of a calculator
    procedure :: setGeometry => TDftbPlus_setGeometry
    !> add an external potential to a calculator
    procedure :: setExternalPotential => TDftbPlus_setExternalPotential
    !> add external charges to a calculator
    procedure :: setExternalCharges => TDftbPlus_setExternalCharges
    !> add reactive external charges to a calculator
    procedure :: setQDepExtPotGen => TDftbPlus_setQDepExtPotGen
    !> obtain the DFTB+ energy
    procedure :: getEnergy => TDftbPlus_getEnergy
    !> obtain the DFTB+ gradients
    procedure :: getGradients => TDftbPlus_getGradients
    !> obtain the DFTB+ stress tensor
    procedure :: getStressTensor => TDftbPlus_getStressTensor
    !> obtain the gradients of the external charges
    procedure :: getExtChargeGradients => TDftbPlus_getExtChargeGradients
    !> get the gross (Mulliken) DFTB+ charges
    procedure :: getGrossCharges => TDftbPlus_getGrossCharges
    !> get the CM5 DFTB+ charges
    procedure :: getCM5Charges => TDftbPlus_getCM5Charges
    !> get electrostatic potential at specified points
    procedure :: getElStatPotential => TDftbPlus_getElStatPotential
    !> Return the number of DFTB+ atoms in the system
    procedure :: nrOfAtoms => TDftbPlus_nrOfAtoms
    !> Return the number of k-points in the DFTB+ calculation (1 if non-repeating)
    procedure :: nrOfKPoints => TDftbPlus_nrOfKPoints
    !> Check that the list of species names has not changed
    procedure :: checkSpeciesNames => TDftbPlus_checkSpeciesNames
    !> Replace species and redefine all quantities that depend on it
    procedure :: setSpeciesAndDependents => TDftbPlus_setSpeciesAndDependents
    !> Initialise electron and nuclear Ehrenfest dynamics
    procedure :: initializeTimeProp => TDftbPlus_initializeTimeProp
    !> Do one propagator step for electrons and, if enabled, nuclei
    procedure :: doOneTdStep => TDftbPlus_doOneTdStep
    !> Set electric field for current propagation step of electrons and nuclei
    procedure :: setTdElectricField => TDftbPlus_setTdElectricField
    !> Set electric field for current propagation step of electrons and nuclei
    procedure :: setTdCoordsAndVelos => TDftbPlus_setTdCoordsAndVelos
    !> Set electric field for current propagation step of electrons and nuclei
    procedure :: getTdForces => TDftbPlus_getTdForces
    !> Check instance of DFTB+ is initialised
    procedure, private :: checkInit => TDftbPlus_checkInit
    !> Return the masses for each atom in the system
    procedure :: getAtomicMasses => TDftbPlus_getAtomicMasses
    !> Return the number of basis functions for each atom in the system
    procedure :: getNOrbitalsOnAtoms => TDftbPlus_getNOrbAtoms
  end type TDftbPlus


#:if not INSTANCE_SAFE_BUILD

  !> Nr. of existing instances (if build is not instance safe)
  integer :: nInstance_ = 0

#:endif


contains


  !> Return the version string for the current DFTB+ build
  subroutine getDftbPlusBuild(version)

    !> Version string for DFTB+
    character(:), allocatable, intent(out) :: version

    version = '${RELEASE}$'

  end subroutine getDftbPlusBuild


  !> Returns the DFTB+ API version
  subroutine getDftbPlusApi(major, minor, patch, instanceSafe)

    !> Major version number
    integer, intent(out) :: major

    !> Minor version number
    integer, intent(out) :: minor

    !> patch level for API
    integer, intent(out) :: patch

    !> Whether API is instance safe
    logical, optional, intent(out) :: instanceSafe

    major = ${APIMAJOR}$
    minor = ${APIMINOR}$
    patch = ${APIPATCH}$
    if (present(instanceSafe)) then
      instanceSafe = instanceSafeBuild
    end if

  end subroutine getDftbPlusApi


  !> Returns the root node of the input, so that it can be further processed
  subroutine TDftbPlusInput_getRootNode(this, root)

    !> Instance.
    class(TDftbPlusInput), intent(in) :: this

    !> Pointer to root node
    type(fnode), pointer, intent(out) :: root

    if (.not. associated(this%hsdTree)) then
      call error("Input has not been created yet!")
    end if
    call getChild(this%hsdTree, rootTag, root)

  end subroutine TDftbPlusInput_getRootNode


  !> Passes the information about the QM region to DFTB+
  subroutine TDftbPlusAtomList_get(instance, nAtom, speciesNames, species)

    !> Input containing the tree representation of the parsed HSD file.
    class(TDftbPlusAtomList), intent(out) :: instance

    !> Number of atoms
    integer, intent(in) :: nAtom

    !> Linked list of chemical symbols of elements (species names), size=nSpecies
    type(TListString), intent(inout) :: speciesNames

    !> Array of species for each atom, size=nAtom
    integer, intent(in) :: species(:)

    integer :: i
    character(3) :: s

    instance%nAtom = nAtom

    call init(instance%speciesNames)
    do i=1,len(speciesNames)
      call get(speciesNames, s, i)
      call append(instance%speciesNames, s)
    end do

    allocate(instance%species(nAtom))
    instance%species(1:nAtom) = species(1:nAtom)

  end subroutine TDftbPlusAtomList_get


  !> Insert the list of atoms into the input data structure
  subroutine TDftbPlusAtomList_addToInpData(instance, inpData)

    !> Input structure of the API
    class(TDftbPlusAtomList), intent(inout) :: instance

    !> Input data structure that will be in turn filled by parsing the HSD tree
    type(TInputData), intent(out), target :: inpData

    type(TGeometry), pointer :: geo

    ! adopted from subroutine readTGeometryGen_help
    geo => inpData%geom

    geo%nAtom = instance%nAtom
    geo%tPeriodic = .false.
    geo%tFracCoord = .false.
    geo%tHelical = .false.

    geo%nSpecies = len(instance%speciesNames)
    allocate(geo%speciesNames(geo%nSpecies))
    call asArray(instance%speciesNames, geo%speciesNames)

    ! Read in sequential and species indices.
    allocate(geo%species(geo%nAtom))
    allocate(geo%coords(3, geo%nAtom))

    geo%species(1:geo%nAtom) = instance%species(1:geo%nAtom)

    if (geo%nSpecies /= maxval(geo%species) .or. minval(geo%species) /= 1) then
      call error("Nr. of species and nr. of specified elements do not match.")
    end if

  end subroutine TDftbPlusAtomList_addToInpData


  !> Initialises a DFTB+ instance
  !>
  !> Note: due to some remaining global variables in the DFTB+ core, only one instance can be
  !> initialised within one process. Therefore, this routine can not be called twice, unless the
  !> TDftbPlus_destruct() has been called in between. Otherwise the subroutine will stop.
  !>
  subroutine TDftbPlus_init(this, outputUnit, mpiComm, devNull)

    !> Instance
    type(TDftbPlus), intent(out) :: this

    !> Unit where to write the output (note: also error messages are written here)
    integer, intent(in), optional :: outputUnit

    !> MPI-communicator to use
    integer, intent(in), optional :: mpiComm

    !> Unit of the null device (you must open the null device and pass its unit number, if you open
    !> multiple TDftbPlus instances within an MPI-process)
    integer, intent(in), optional :: devNull

    integer :: stdOut

    #:if not INSTANCE_SAFE_BUILD
      if (nInstance_ /= 0) then
        call error("This build does not support multiple DFTB+ instances")
      end if
      nInstance_ = 1
    #:endif

    if (present(mpiComm) .and. .not. withMpi) then
      call error("MPI Communicator supplied to initialise serial DFTB+ instance")
    end if

    if (present(outputUnit)) then
      stdOut = outputUnit
    else
      stdOut = output_unit
    end if

    call initGlobalEnv(outputUnit=outputUnit, mpiComm=mpiComm, devNull=devNull)
    allocate(this%env)
    allocate(this%main)
    call TEnvironment_init(this%env)
    this%env%tAPICalculation = .true.
    this%tInit = .true.

  end subroutine TDftbPlus_init


  !> Destroys a DFTB+ instance
  subroutine TDftbPlus_destruct(this)

    !> Instance
    type(TDftbPlus), intent(inout) :: this

    call this%checkInit()

    call this%main%destructProgramVariables()
    call this%env%destruct()
    deallocate(this%main, this%env)
    call destructGlobalEnv()
    this%tInit = .false.

    #:if not INSTANCE_SAFE_BUILD
      nInstance_ = 0
    #:endif

  end subroutine TDftbPlus_destruct


  !> Fills up the input by parsing an HSD file
  subroutine TDftbPlus_getInputFromFile(this, fileName, input)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Name of the file to parse
    character(len=*), intent(in) :: fileName

    !> Input containing the tree representation of the parsed HSD file.
    type(TDftbPlusInput), intent(out) :: input

    call this%checkInit()

    call readHsdFile(fileName, input%hsdTree)

  end subroutine TDftbPlus_getInputFromFile


  !> Creates an input with no entries.
  subroutine TDftbPlus_getEmptyInput(this, input)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Instance.
    type(TDftbPlusInput), intent(out) :: input

    type(fnode), pointer :: root, dummy

    call this%checkInit()

    input%hsdTree => createDocumentNode()
    root => createElement(rootTag)
    dummy => appendChild(input%hsdTree, root)

  end subroutine TDftbPlus_getEmptyInput


  !> Sets up the calculator using a given input.
  subroutine TDftbPlus_setupCalculator(this, input)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Representation of the DFTB+ input.
    type(TDftbPlusInput), intent(inout) :: input

    type(TParserFlags) :: parserFlags
    type(TInputData) :: inpData

    call this%checkInit()

    call parseHsdTree(input%hsdTree, inpData, parserFlags)
    call doPostParseJobs(input%hsdTree, parserFlags)
    call this%main%initProgramVariables(inpData, this%env)

  end subroutine TDftbPlus_setupCalculator


  !> Sets the geometry in the calculator.
  subroutine TDftbPlus_setGeometry(this, coords, latVecs, origin)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Atomic coordinates in Bohr units. Shape: (3, nAtom).
    real(dp), intent(in) :: coords(:,:)

    !> Lattice vectors in Bohr units, stored column-wise. Shape: (3, 3).
    real(dp), intent(in), optional :: latVecs(:,:)

    !> Coordinate origin in Bohr units. Shape: (3).
    real(dp), intent(in), optional :: origin(:)

    call this%checkInit()

    call setGeometry(this%env, this%main, coords, latVecs, origin)

  end subroutine TDftbPlus_setGeometry


  !> Sets an external potential.
  subroutine TDftbPlus_setExternalPotential(this, atomPot, potGrad)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Potential acting on each atom. Shape: (nAtom)
    real(dp), intent(in), optional :: atomPot(:)

    !> Gradient of the potential  on each atom. Shape: (3, nAtom)
    real(dp), intent(in), optional :: potGrad(:,:)

    call this%checkInit()

    if (allocated(this%main%solvation)) then
      if (this%main%solvation%isEFieldModified()) then
        call error("External fields currently unsupported for this solvent model")
      end if
    end if

    call setExternalPotential(this%main, atomPot=atomPot, potGrad=potGrad)

  end subroutine TDftbPlus_setExternalPotential


  !> Sets external point charges.
  subroutine TDftbPlus_setExternalCharges(this, chargeCoords, chargeQs, blurWidths)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Coordinate of the external charges
    real(dp), intent(in) :: chargeCoords(:,:)

    !> Charges of the external point charges (sign convention: electron is negative)
    real(dp), intent(in) :: chargeQs(:)

    !> Widths of the Gaussian for each charge used for blurring (0.0 = no blurring)
    real(dp), intent(in), optional :: blurWidths(:)

    call this%checkInit()

    if (allocated(this%main%solvation)) then
      if (this%main%solvation%isEFieldModified()) then
        call error("External fields currently unsupported for this solvent model")
      end if
    end if

    call setExternalCharges(this%main, chargeCoords, chargeQs, blurWidths)

  end subroutine TDftbPlus_setExternalCharges


  !> Sets the generator for the population dependant external potential.
  subroutine TDftbPlus_setQDepExtPotGen(this, extPotGen)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Population dependant external potential generator
    class(TQDepExtPotGen), intent(in) :: extPotGen

    type(TQDepExtPotGenWrapper) :: extPotGenWrapper
    type(TQDepExtPotProxy) :: extPotProxy

    call this%checkInit()

    if (allocated(this%main%solvation)) then
      if (this%main%solvation%isEFieldModified()) then
        call error("External fields currently unsupported for this solvent model")
      end if
    end if

    allocate(extPotGenWrapper%instance, source=extPotGen)
    call TQDepExtPotProxy_init(extPotProxy, [extPotGenWrapper])
    call setQDepExtPotProxy(this%main, extPotProxy)

  end subroutine TDftbPlus_setQDepExtPotGen


  !> Return the energy of the current system.
  subroutine TDftbPlus_getEnergy(this, merminEnergy)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Mermin free energy.
    real(dp), intent(out) :: merminEnergy

    call this%checkInit()

    call getEnergy(this%env, this%main, merminEnergy)

  end subroutine TDftbPlus_getEnergy


  !> Returns the gradient on the atoms in the system.
  subroutine TDftbPlus_getGradients(this, gradients)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Gradients on the atoms.
    real(dp), intent(out) :: gradients(:,:)

    call this%checkInit()

    call getGradients(this%env, this%main, gradients)

  end subroutine TDftbPlus_getGradients


  !> Returns the stress tensor of the periodic system.
  subroutine TDftbPlus_getStressTensor(this, stresstensor)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Gradients on the atoms.
    real(dp), intent(out) :: stresstensor(:,:)

    call this%checkInit()

    call getStressTensor(this%env, this%main, stresstensor)

  end subroutine TDftbPlus_getStressTensor


  !> Returns the gradients on the external charges.
  !>
  !> This function may only be called if TDftbPlus_setExternalCharges was called before it
  !>
  subroutine TDftbPlus_getExtChargeGradients(this, gradients)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Gradients acting on the external charges.
    real(dp), intent(out) :: gradients(:,:)

    call this%checkInit()

    call getExtChargeGradients(this%main, gradients)

  end subroutine TDftbPlus_getExtChargeGradients


  !> Returns the gross charges of each atom
  subroutine TDftbPlus_getGrossCharges(this, atomCharges)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Atomic gross charges.
    real(dp), intent(out) :: atomCharges(:)

    call this%checkInit()

    call getGrossCharges(this%env, this%main, atomCharges)

  end subroutine TDftbPlus_getGrossCharges


  !> Returns the CM5 charges of each atom
  subroutine TDftbPlus_getCM5Charges(this, atomCharges)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Atomic gross charges.
    real(dp), intent(out) :: atomCharges(:)

    call this%checkInit()

    call getCM5Charges(this%env, this%main, atomCharges)

  end subroutine TDftbPlus_getCM5Charges


  !> Returns electrostatic potential at specified points
  subroutine TDftbPlus_getElStatPotential(this, pot, locations)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Resulting potentials
    real(dp), intent(out) :: pot(:)

    !> Sites at which to calculate potential
    real(dp), intent(in) :: locations(:,:)

    call this%checkInit()

    call getElStatPotential(this%env, this%main, pot, locations)

  end subroutine TDftbPlus_getElStatPotential


  !> Returns the nr. of atoms in the system.
  function TDftbPlus_nrOfAtoms(this) result(nAtom)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Nr. of atoms
    integer :: nAtom

    call this%checkInit()

    nAtom = nrOfAtoms(this%main)

  end function TDftbPlus_nrOfAtoms


  !> Returns the nr. of k-points describing the system.
  function TDftbPlus_nrOfKPoints(this) result(nKpts)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Nr. of k-points
    integer :: nKpts

    call this%checkInit()

    nKpts = nrOfKPoints(this%main)

  end function TDftbPlus_nrOfKPoints


  !> Returns the atomic masses for each atom in the system.
  subroutine TDftbPlus_getAtomicMasses(this, mass)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Masses for each species of the system
    real(dp), intent(out) :: mass(:)

    call this%checkInit()

    call getAtomicMasses(this%main, mass)

  end subroutine TDftbPlus_getAtomicMasses


  !> Returns the number of orbitals for each atom in the system
  subroutine TDftbPlus_getNOrbAtoms(this, nOrbs)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Number of basis functions associated with each atom
    integer, intent(out) :: nOrbs(:)

    nOrbs(:) = this%main%orb%nOrbAtom

  end subroutine TDftbPlus_getNOrbAtoms


  !> Checks whether the type is already initialized and stops the code if not.
  subroutine TDftbPlus_checkInit(this)

    !> Instance.
    class(TDftbPlus), intent(in) :: this

    if (.not. this%tInit) then
      call error("Received uninitialized TDftbPlus instance")
    end if

  end subroutine TDftbPlus_checkInit


  !> Reads out the atomic angular momenta from an SK-file
  !>
  !> NOTE: This only works with handcrafted (non-standard) SK-files, where the nr. of shells
  !>   has been added as 3rd entry to the first line of the homo-nuclear SK-files.
  !>
  function getMaxAngFromSlakoFile(slakoFile) result(maxAng)

    !> Instance.
    character(len=*), intent(in) :: slakoFile

    !> Maximal angular momentum found in the file
    integer :: maxAng

    real(dp) :: dr
    integer :: nGridPoints, nShells
    integer :: fd

    open(newunit=fd, file=slakoFile, status="old", action="read")
    read(fd, *) dr, nGridPoints, nShells
    close(fd)
    maxAng = nShells - 1

  end function getMaxAngFromSlakoFile


  !> Converts atom types to species
  subroutine convertAtomTypesToSpecies(typeNumbers, species, speciesNames, typeNames)

    !> Type number of each atom is the system. It can be arbitrary number (e.g. atomic number)
    integer, intent(in) :: typeNumbers(:)

    !> Species index for each atom (1 for the first atom type found, 2 for the second, etc.)
    integer, allocatable, intent(out) :: species(:)

    !> Names of each species, usually X1, X2 unless typeNames have been specified.
    character(len=*), allocatable, intent(out) :: speciesNames(:)

    !> Array of type names, indexed by the type numbers.
    character(len=*), intent(in), optional :: typeNames(:)

    integer, allocatable :: uniqueTypeNumbers(:)
    integer :: nAtom, nSpecies
    integer :: iAt, iSp, curType

    nAtom = size(typeNumbers)

    allocate(uniqueTypeNumbers(nAtom))
    nSpecies = 0
    do iAt = 1, nAtom
      curType = typeNumbers(iAt)
      if (.not. any(uniqueTypeNumbers(1:nSpecies) == curType)) then
        nSpecies = nSpecies + 1
        uniqueTypeNumbers(nSpecies) = curType
      end if
    end do

    allocate(species(nAtom))
    do iSp = 1, nSpecies
      where (typeNumbers == uniqueTypeNumbers(iSp))
        species = iSp
      end where
    end do

    allocate(speciesNames(nSpecies))
    do iSp = 1, nSpecies
      if (present(typeNames)) then
        speciesNames(iSp) = typeNames(uniqueTypeNumbers(iSp))
      else
        write(speciesNames(iSp), "(A,I0)") "X", iSp
      end if
    end do

  end subroutine convertAtomTypesToSpecies


  !> Check whether speciesNames has changed between calls to DFTB+
  subroutine TDftbPlus_checkSpeciesNames(this, inputSpeciesNames)
    class(TDftbPlus),  intent(inout) :: this
    character(len=*), intent(in) :: inputSpeciesNames(:)

    logical :: tSpeciesNameChanged

    call this%checkInit()

    tSpeciesNameChanged = checkSpeciesNames(this%env, this%main, inputSpeciesNames)

    if(tSpeciesNameChanged)then
      call error('speciesNames has changed between calls to DFTB+. This will cause erroneous&
          & results.' // newline // 'Instead call destruct and then fully re-initialize.')
    else
       continue
    endif

  end subroutine TDftbPlus_checkSpeciesNames


  !> Set species and all variables/data dependent on it
  subroutine TDftbPlus_setSpeciesAndDependents(this, inputSpeciesNames, inputSpecies)
    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Type of each atom (nAllAtom)
    integer, intent(in) :: inputSpecies(:)

    !> Labels of atomic species (nSpecies)
    character(len=*), intent(in) :: inputSpeciesNames(:)

    call this%checkInit()
    call this%checkSpeciesNames(inputSpeciesNames)
    call updateDataDependentOnSpeciesOrdering(this%env, this%main, inputSpecies)

  end subroutine TDftbPlus_setSpeciesAndDependents


  !> Initialise propagators for electron and nuclei dynamics
  subroutine TDftbPlus_initializeTimeProp(this, dt, tdFieldThroughAPI, tdCoordsAndVelosThroughAPI)
    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> time step
    real(dp), intent(in) :: dt

    !> field will be provided through the API?
    logical, intent(in) :: tdFieldThroughAPI

    !> coords and velocities will be provided at each step through the API?
    logical, intent(in) :: tdCoordsAndVelosThroughAPI

    call initializeTimeProp(this%env, this%main, dt, tdFieldThroughAPI, tdCoordsAndVelosThroughAPI)

  end subroutine TDftbPlus_initializeTimeProp


  !> Propagate one time step for electron and nuclei dynamics
  subroutine TDftbPlus_doOneTdStep(this, iStep, dipole, energy, atomNetCharges,&
      & coord, force, occ)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> present step of dynamics
    integer, intent(in) :: iStep

    !> Dipole moment
    real(dp), optional, intent(out) :: dipole(:,:)

    !> data type for energy components and total
    real(dp), optional, intent(out) :: energy

    !> Negative gross charge
    real(dp), optional, intent(out) :: atomNetCharges(:,:)

    !> atomic coordinates
    real(dp), optional, intent(out) :: coord(:,:)

    !> forces (3, nAtom)
    real(dp), optional, intent(out) :: force(:,:)

    !> molecular orbital projected populations
    real(dp), optional, intent(out) :: occ(:)

    call doOneTdStep(this%env, this%main, iStep, dipole=dipole, energy=energy,&
        & atomNetCharges=atomNetCharges, coordOut = coord, force=force, occ=occ)

  end subroutine TDftbPlus_doOneTdStep


  !> sets electric field for td propagation
  subroutine TDftbPlus_setTdElectricField(this, field)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> electric field components
    real(dp), intent(in) :: field(3)

    if (allocated(this%main%solvation)) then
      if (this%main%solvation%isEFieldModified()) then
        call error("External fields currently unsupported for this solvent model")
      end if
    end if

    call setTdElectricField(this%main, field)

  end subroutine TDftbPlus_setTdElectricField


  !> Set atomic coordinates and velocities for MD
  subroutine TDftbPlus_setTdCoordsAndVelos(this, coords, velos)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> coordinates
    real(dp), intent(in) :: coords(3, this%main%nAtom)

    !> velocities
    real(dp), intent(in) :: velos(3, this%main%nAtom)

    call setTdCoordsAndVelos(this%main, coords, velos)

  end subroutine TDftbPlus_setTdCoordsAndVelos


  !> Returns forces from time dependent propagation
  subroutine TDftbPlus_getTdForces(this, forces)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> forces (3, nAtom)
    real(dp), intent(out) :: forces(:,:)

    call getTdForces(this%main, forces)

  end subroutine TDftbPlus_getTdForces

end module dftbp_mmapi
