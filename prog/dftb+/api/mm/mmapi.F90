!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Proviedes DFTB+ API for MM-type high level access
module dftbp_mmapi
  use iso_fortran_env, only : output_unit
  use dftbp_accuracy
  use dftbp_globalenv
  use dftbp_environment
  use dftbp_message, only: error
  use dftbp_mainapi
  use dftbp_parser
  use dftbp_hsdutils
  use dftbp_hsdhelpers, only : doPostParseJobs
  use dftbp_inputdata
  use dftbp_xmlf90
  use dftbp_qdepextpotgen, only : TQDepExtPotGen, TQDepExtPotGenWrapper
  use dftbp_qdepextpotproxy, only : TQDepExtPotProxy, TQDepExtPotProxy_init
  use dftbp_charmanip, only : newline
  implicit none
  private

  public :: TDftbPlus, getDftbPlusBuild, getDftbPlusApi
  public :: TDftbPlus_init, TDftbPlus_destruct
  public :: TDftbPlusInput
  public :: TQDepExtPotGen
  public :: getMaxAngFromSlakoFile, convertAtomTypesToSpecies


  !> Number of DFTB+ calculation instances
  integer :: nDftbPlusCalc = 0


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
    type(TEnvironment) :: env
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
    !> Return the number of DFTB+ atoms in the system
    procedure :: nrOfAtoms => TDftbPlus_nrOfAtoms
    !> Check that the list of species names has not changed
    procedure :: checkSpeciesNames => TDftbPlus_checkSpeciesNames
    !> Replace species and redefine all quantities that depend on it
    procedure :: setSpeciesAndDependents => TDftbPlus_setSpeciesAndDependents
    !> Check instance of DFTB+ is initialised
    procedure, private :: checkInit => TDftbPlus_checkInit
  end type TDftbPlus



contains


  !> Return the version string for the current DFTB+ build
  subroutine getDftbPlusBuild(version)

    !> Version string for DFTB+
    character(:), allocatable, intent(out) :: version

    version = '${RELEASE}$'

  end subroutine getDftbPlusBuild


  !> Returns the DFTB+ API version
  subroutine getDftbPlusApi(major, minor, patch)

    !> Major version number
    integer, intent(out) :: major

    !> Minor version number
    integer, intent(out) :: minor

    !> patch level for API
    integer, intent(out) :: patch

    major = ${APIMAJOR}$
    minor = ${APIMINOR}$
    patch = ${APIPATCH}$

  end subroutine getDftbPlusApi


  !> Returns the root node of the input, so that it can be further processed
  subroutine TDftbPlusInput_getRootNode(this, root)

    !> Instance.
    class(TDftbPlusInput), intent(in) :: this

    !> Pointer to root node
    type(fnode), pointer, intent(out) :: root

    if (.not. associated(this%hsdTree)) then
      print *, 'ERROR: input has not been created yet!'
      stop
    end if
    call getChild(this%hsdTree, rootTag, root)

  end subroutine TDftbPlusInput_getRootNode


  !> Initalises a DFTB+ instance
  !>
  !> Note: due to some remaining global variables in the DFTB+ core, only one instance can be
  !> initialised within one process. Therefore, this routine can not be called twice, unless the
  !> TDftbPlus_destruct() has been called in between. Otherwise the subroutine will stop.
  !>
  subroutine TDftbPlus_init(this, outputUnit, mpiComm)

    !> Instance
    type(TDftbPlus), intent(out) :: this

    !> Unit where to write the output (note: also error messages are written here)
    integer, intent(in), optional :: outputUnit

    !> MPI-communicator to use
    integer, intent(in), optional :: mpiComm

    integer :: stdOut

    if (present(outputUnit)) then
      stdOut = outputUnit
    else
      stdOut = output_unit
    end if

    if (nDftbPlusCalc /= 0) then
      write(stdOut, "(A)") "Error: only one DFTB+ instance is currently allowed"
      stop
    end if
    nDftbPlusCalc = 1

    call initGlobalEnv(outputUnit=outputUnit, mpiComm=mpiComm)
    call TEnvironment_init(this%env)
    this%env%tAPICalculation = .true.
    this%tInit = .true.

  end subroutine TDftbPlus_init


  !> Destroys a DFTB+ instance
  subroutine TDftbPlus_destruct(this)

    !> Instance
    type(TDftbPlus), intent(inout) :: this

    call this%checkInit()

    call destructProgramVariables()
    call this%env%destruct()
    call destructGlobalEnv()
    nDftbPlusCalc = 0
    this%tInit = .false.

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
    call initProgramVariables(inpData, this%env)

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

    call setGeometry(this%env, coords, latVecs, origin)

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

    call setExternalPotential(atomPot=atomPot, potGrad=potGrad)

  end subroutine TDftbPlus_setExternalPotential


  !> Sets external point charges.
  subroutine TDftbPlus_setExternalCharges(this, chargeCoords, chargeQs, blurWidths)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Coordiante of the external charges
    real(dp), intent(in) :: chargeCoords(:,:)

    !> Charges of the external point charges (sign convention: electron is negative)
    real(dp), intent(in) :: chargeQs(:)

    !> Widths of the Gaussian for each charge used for blurring (0.0 = no blurring)
    real(dp), intent(in), optional :: blurWidths(:)

    call this%checkInit()

    call setExternalCharges(chargeCoords, chargeQs, blurWidths)

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

    allocate(extPotGenWrapper%instance, source=extPotGen)
    call TQDepExtPotProxy_init(extPotProxy, [extPotGenWrapper])
    call setQDepExtPotProxy(extPotProxy)

  end subroutine TDftbPlus_setQDepExtPotGen


  !> Return the energy of the current system.
  subroutine TDftbPlus_getEnergy(this, merminEnergy)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Mermin free energy.
    real(dp), intent(out) :: merminEnergy

    call this%checkInit()

    call getEnergy(this%env, merminEnergy)

  end subroutine TDftbPlus_getEnergy


  !> Returns the gradient on the atoms in the system.
  subroutine TDftbPlus_getGradients(this, gradients)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Gradients on the atoms.
    real(dp), intent(out) :: gradients(:,:)

    call this%checkInit()

    call getGradients(this%env, gradients)

  end subroutine TDftbPlus_getGradients


     !> Returns the stress tensor of the periodic system.
  subroutine TDftbPlus_getStressTensor(this, stresstensor)

    !> Instance.
    class(TDftbPlus), intent(inout) :: this

    !> Gradients on the atoms.
    real(dp), intent(out) :: stresstensor(:,:)

    call this%checkInit()

    call getStressTensor(this%env, stresstensor)

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

    call getExtChargeGradients(gradients)

  end subroutine TDftbPlus_getExtChargeGradients


  !> Returns the gross charges of each atom
  subroutine TDftbPlus_getGrossCharges(this, atomCharges)

    !> Instance
    class(TDftbPlus), intent(inout) :: this

    !> Atomic gross charges.
    real(dp), intent(out) :: atomCharges(:)

    call this%checkInit()

    call getGrossCharges(this%env, atomCharges)

  end subroutine TDftbPlus_getGrossCharges


  !> Returns the nr. of atoms in the system.
  function TDftbPlus_nrOfAtoms(this) result(nAtom)

    !> Instance
    class(TDftbPlus), intent(in) :: this

    !> Nr. of atoms
    integer :: nAtom

    call this%checkInit()

    nAtom = nrOfAtoms()

  end function TDftbPlus_nrOfAtoms


  !> Checks whether the type is already initialized and stops the code if not.
  subroutine TDftbPlus_checkInit(this)

    !> Instance.
    class(TDftbPlus), intent(in) :: this

    if (.not. this%tInit) then
      write(stdOut, "(A)") "Error: Received uninitialized TDftbPlus instance"
      stop
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

    tSpeciesNameChanged = checkSpeciesNames(this%env, inputSpeciesNames)

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
    call updateDataDependentOnSpeciesOrdering(this%env, inputSpecies)

  end subroutine TDftbPlus_setSpeciesAndDependents

end module dftbp_mmapi
