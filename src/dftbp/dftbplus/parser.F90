!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Fills the derived type with the input parameters from an HSD or an XML file.
module dftbp_dftbplus_parser
  use dftbp_common_accuracy, only : dp, sc, lc, mc, minTemp, distFudge, distFudgeOld
  use dftbp_common_constants, only : pi, boltzmann, Bohr__AA, maxL, shellNames, symbolToNumber
  use dftbp_common_file, only : fileAccessValues, openFile, closeFile, TFileDescr
  use dftbp_common_filesystem, only : findFile, getParamSearchPath
  use dftbp_common_globalenv, only : stdout, withMpi, withScalapack, abortProgram
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_release, only : TVersionMap
  use dftbp_common_status, only : TStatus
  use dftbp_common_unitconversion, only : lengthUnits, energyUnits, forceUnits, pressureUnits,&
      & timeUnits, EFieldUnits, freqUnits, massUnits, VelocityUnits, dipoleUnits, chargeUnits,&
      & volumeUnits, angularUnits
  use dftbp_dftb_elecconstraints, only : readElecConstraintInput
  use dftbp_dftb_coordnumber, only : TCNInput, getElectronegativity, getCovalentRadius, cnType
  use dftbp_dftb_dftbplusu, only : plusUFunctionals
  use dftbp_dftb_dftd4param, only : getEeqChi, getEeqGam, getEeqKcn, getEeqRad
  use dftbp_dftb_dispersions, only : TDispersionInp, TDispSlaKirkInp, TDispUffInp,&
      & TSimpleDftD3Input, TDispDftD4Inp, getUffValues
  use dftbp_dftb_encharges, only : TEeqInput
  use dftbp_dftb_etemp, only : fillingTypes
  use dftbp_dftb_halogenx, only : halogenXSpecies1, halogenXSpecies2
  use dftbp_dftb_periodic, only : TNeighbourList, TNeighbourlist_init, getSuperSampling, &
      & getCellTranslations, updateNeighbourList
  use dftbp_dftb_rangeseparated, only : TRangeSepSKTag, rangeSepTypes, rangeSepFunc
  use dftbp_dftb_repulsive_chimesrep, only : TChimesRepInp
  use dftbp_dftb_repulsive_polyrep, only : TPolyRepInp, TPolyRep
  use dftbp_dftb_repulsive_splinerep, only : TSplineRepInp, TSplineRep
  use dftbp_dftb_slakocont, only : init, addTable
  use dftbp_dftb_slakoeqgrid, only : skEqGridNew, skEqGridOld, TSlakoEqGrid, init
  use dftbp_dftbplus_forcetypes, only : forceTypes
  use dftbp_dftbplus_input_fileaccess, only : readBinaryAccessTypes
  use dftbp_dftbplus_inputconversion, only : transformpdosregioninfo
  use dftbp_dftbplus_inputdata, only :TInputData, TControl, TSlater, TBlacsOpts, TRangeSepInp
  use dftbp_dftbplus_oldcompat, only : convertOldHSD
  use dftbp_dftbplus_specieslist, only : readSpeciesList
  use dftbp_elecsolvers_elecsolvers, only : electronicSolverTypes
  use dftbp_extlibs_arpack, only : withArpack
  use dftbp_extlibs_elsiiface, only : withELSI, withPEXSI
  use dftbp_extlibs_plumed, only : withPlumed
  use dftbp_extlibs_poisson, only : withPoisson, TPoissonInfo, TPoissonStructure
  use dftbp_extlibs_sdftd3, only : TSDFTD3Input, dampingFunction
  use dftbp_extlibs_tblite, only : tbliteMethod
  use dftbp_extlibs_xmlf90, only : fnode, removeChild, string, char, textNodeName, fnodeList,&
      & getLength, getNodeName, getItem1, destroyNodeList, destroyNode, assignment(=)
  use dftbp_geoopt_geoopt, only : geoOptTypes
  use dftbp_dftbplus_input_geoopt, only : readGeoOptInput
  use dftbp_io_charmanip, only : i2c, newline, tolower, unquote
  use dftbp_io_hsdparser, only : getNodeHSdName, parseHsd
  use dftbp_io_hsdutils, only : detailedError, detailedWarning, getChild, getChildValue,&
      & getChildren, getSelectedAtomIndices, setChild, setChildValue
  use dftbp_io_hsdutils2, only : convertUnitHsd, getNodeName2, setUnprocessed, splitModifier,&
      & renameChildren
  use dftbp_io_message, only : error, warning
  use dftbp_io_xmlutils, only : removeChildNodes
  use dftbp_math_lapackroutines, only : matinv
  use dftbp_math_simplealgebra, only: cross3, determinant33
  use dftbp_md_tempprofile, only : identifyTempProfile
  use dftbp_md_xlbomd, only : TXlbomdInp
  use dftbp_mixer_mixer, only : mixerTypes
  use dftbp_dftb_nonscc, only : diffTypes
  use dftbp_reks_reks, only : reksTypes
  use dftbp_solvation_solvparser, only : readSolvation, readCM5
  use dftbp_timedep_linresptypes, only : linRespSolverTypes
  use dftbp_timedep_timeprop, only : TElecDynamicsInp, pertTypes, tdSpinTypes, envTypes
  use dftbp_type_commontypes, only : TOrbitals
  use dftbp_type_linkedlist, only : TListCharLc, TListInt, TListIntR1, TListReal, TListRealR1,&
      & TListRealR2, TListString, init, destruct, append, get, len, asArray, asVector, intoArray
  use dftbp_type_oldskdata, only : TOldSKData, readFromFile, inquireRangeSepTag
  use dftbp_type_orbitals, only : getShellnames
  use dftbp_type_typegeometry, only : TGeometry, reduce, setLattice
  use dftbp_type_typegeometryhsd, only : readTGeometryGen, readTGeometryHsd, readTGeometryXyz,&
      & readTGeometryVasp, readTGeometryLammps
  use dftbp_type_wrappedintr, only : TWrappedInt1
#:if WITH_MBD
  use dftbp_dftb_dispmbd, only :TDispMbdInp
#:endif
#:if WITH_SOCKETS
  use dftbp_io_ipisocket, only : IPI_PROTOCOLS
#:endif
#:if WITH_TRANSPORT
  use dftbp_transport_negfvars, only : TTransPar, TNEGFGreenDensInfo, TNEGFTunDos, TElPh,&
      & ContactInfo
#:endif
  implicit none

  private
  public :: parserVersion, rootTag
  public :: TParserFlags
  public :: readHsdFile, parseHsdTree


  !> Tag at the head of the input document tree
  character(len=*), parameter :: rootTag = "dftbplusinput"


  !> Container type for parser related flags.
  type TParserFlags

    !> stop after parsing?
    logical :: tStop

    !> Continue despite unprocessed nodes
    logical :: tIgnoreUnprocessed

    !> HSD output?
    logical :: tWriteHSD
  end type TParserFlags

  !> Actual input version <-> parser version maps (must be updated at every public release)
  type(TVersionMap), parameter :: versionMaps(*) = [&
      & TVersionMap("23.2", 14), TVersionMap("23.1", 13), TVersionMap("22.2", 12),&
      & TVersionMap("22.1", 11), TVersionMap("21.2", 10), TVersionMap("21.1", 9),&
      & TVersionMap("20.2", 9), TVersionMap("20.1", 8), TVersionMap("19.1", 7),&
      & TVersionMap("18.2", 6), TVersionMap("18.1", 5), TVersionMap("17.1", 5)]

  !> Version of the oldest parser for which compatibility is still maintained
  integer, parameter :: minVersion = 1

  !> Version of the current parser (as latest version)
  integer, parameter :: parserVersion = maxval(versionMaps(:)%parserVersion)


contains

  !> Reads the HSD input from a file
  subroutine readHsdFile(hsdFile, hsdTree, errStatus)

    !> Name of the input file
    character(*), intent(in) :: hsdFile

    !> Data tree representation of the input
    type(fnode), pointer :: hsdTree

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    call parseHSD(rootTag, hsdFile, hsdTree, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readHsdFile


  !> Parse input from an HSD/XML file
  subroutine parseHsdTree(hsdTree, input, parserFlags, errStatus)

    !> Tree representation of the input
    type(fnode), pointer :: hsdTree

    !> Returns initialised input variables on exit
    type(TInputData), intent(out) :: input

    !> Special block containings parser related settings
    type(TParserFlags), intent(out) :: parserFlags

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(TOrbitals) :: orb
    type(fnode), pointer :: root, tmp, driverNode, hamNode, analysisNode, child, dummy
    logical :: tReadAnalysis
    integer, allocatable :: implicitParserVersion

    write(stdout, '(A,1X,I0,/)') 'Parser version:', parserVersion
    write(stdout, "(A)") repeat("-", 80)

    call getChild(hsdTree, rootTag, root, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call handleInputVersion(root, implicitParserVersion, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(root, "ParserOptions", dummy, errStatus, "", child=child, list=.true.,&
        & allowEmptyValue=.true., dummyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call readParserOptions(child, root, parserFlags, errStatus, implicitParserVersion)
    @:PROPAGATE_ERROR(errStatus)

    ! Read the geometry unless the list of atoms has been provided through the API
    if (.not. allocated(input%geom%coords)) then
      call getChild(root, "Geometry", tmp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call readGeometry(tmp, input, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! Hamiltonian settings that need to know settings from the REKS block
    call getChildValue(root, "Reks", dummy, errStatus, "None", child=child)
    @:PROPAGATE_ERROR(errStatus)
    call readReks(child, dummy, input%ctrl, input%geom, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChild(root, "Transport", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)

  #:if WITH_TRANSPORT

    ! Read in transport and modify geometry if it is only a contact calculation
    if (associated(child)) then
      call readTransportGeometry(child, input%geom, input%transpar, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else
      input%transpar%ncont=0
      allocate(input%transpar%contacts(0))
      ! set range of atoms in the 'device', as there are no contacts
      input%transpar%idxdevice(1) = 1
      input%transpar%idxdevice(2) = input%geom%nAtom
    end if

    call getChild(root, "Dephasing", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      call detailedError(child, "Be patient... Dephasing feature will be available soon!",&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
      !call readDephasing(child, input%slako%orb, input%geom, input%transpar, input%ginfo%tundos)
    end if

    ! electronic Hamiltonian
    call getChildValue(root, "Hamiltonian", hamNode, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call readHamiltonian(hamNode, input%ctrl, input%geom, input%slako, input%transpar,&
        & input%ginfo%greendens, input%poisson, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  #:else

    if (associated(child)) then
      call detailedError(child, "Program has been compiled without transport enabled", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! electronic Hamiltonian
    call getChildValue(root, "Hamiltonian", hamNode, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call readHamiltonian(hamNode, input%ctrl, input%geom, input%slako, input%poisson, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  #:endif

    call getChildValue(root, "Driver", driverNode, errStatus, "", child=child,&
        & allowEmptyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
  #:if WITH_TRANSPORT
    call readDriver(driverNode, child, input%geom, input%ctrl, input%transpar, errStatus)
    @:PROPAGATE_ERROR(errStatus)
  #:else
    call readDriver(driverNode, child, input%geom, input%ctrl, errStatus)
    @:PROPAGATE_ERROR(errStatus)
  #:endif

    tReadAnalysis = .true.
    call getChild(root, "ElectronDynamics", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      allocate(input%ctrl%elecDynInp)
      call readElecDynamics(child, input%ctrl%elecDynInp, input%geom, input%ctrl%masses, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (input%ctrl%elecDynInp%tReadRestart .and. .not.input%ctrl%elecDynInp%tPopulations) then
        tReadAnalysis = .false.
      end if

    end if

    if (tReadAnalysis) then
      ! Analysis of properties
      call getChildValue(root, "Analysis", dummy, errStatus, "", child=analysisNode, list=.true.,&
          & allowEmptyValue=.true., dummyValue=.true.)
      @:PROPAGATE_ERROR(errStatus)

    #:if WITH_TRANSPORT
      call readAnalysis(analysisNode, input%ctrl, input%geom, input%slako%orb, input%transpar,&
          & input%ginfo%tundos, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call finalizeNegf(input, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:else
      call readAnalysis(analysisNode, input%ctrl, input%geom, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:endif

    end if

    call getChildValue(root, "ExcitedState", dummy, errStatus, "", child=child, list=.true.,&
        & allowEmptyValue=.true., dummyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call readExcited(child, input%geom, input%ctrl, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Hamiltonian settings that need to know about settings from the blocks above
    call readLaterHamiltonian(hamNode, input%ctrl, driverNode, input%geom, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(root, "Options", dummy, errStatus, "", child=child, list=.true.,&
        & allowEmptyValue=.true., dummyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call readOptions(child, input%ctrl, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! W values if needed by Hamiltonian or excited state calculation
    if (allocated(input%ctrl%tbliteInp)) then
      call input%ctrl%tbliteInp%setupOrbitals(input%geom%species, orb)
      call readSpinConstants(hamNode, input%geom, orb, input%ctrl, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else
      call readSpinConstants(hamNode, input%geom, input%slako%orb, input%ctrl, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! analysis settings that need to know settings from the options block
    if (tReadAnalysis) then
      call readLaterAnalysis(analysisNode, input%ctrl, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! read parallel calculation settings
    call readParallel(root, input, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! input data strucutre has been initialised
    input%tInitialized = .true.

  end subroutine parseHsdTree


  !> Converts input version to parser version and removes InputVersion node if present.
  subroutine handleInputVersion(root, implicitParserVersion, errStatus)

    !> Root eventually containing InputVersion
    type(fnode), pointer, intent(in) :: root

    !> Parser version corresponding to input version, or unallocated if none has been found
    integer, allocatable, intent(out) :: implicitParserVersion

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child, dummy
    type(string) :: versionString

    call getChild(root, "InputVersion", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      call getChildValue(child, "", versionString, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      allocate(implicitParserVersion)
      call parserVersionFromInputVersion(unquote(char(versionString)), child,&
          & implicitParserVersion, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      dummy => removeChild(root, child)
      call destroyNode(dummy)
    end if

  end subroutine handleInputVersion


  !> Read in parser options (options not passed to the main code)
  subroutine readParserOptions(node, root, flags, errStatus, implicitVersion)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Root of the entire tree (in case it needs to be converted, for example because of
    !> compatibility options)
    type(fnode), pointer :: root

    !> Contains parser flags on exit.
    type(TParserFlags), intent(out) :: flags

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> Parser version implied by version number
    integer, intent(in), optional :: implicitVersion

    integer :: inputVersion
    type(fnode), pointer :: child

    if (present(implicitVersion)) then
      call getChild(node, "ParserVersion", child, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child)) then
        call getChildValue(child, "", inputVersion, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        if (inputVersion /= implicitVersion) then
          call detailedError(child, "Parser version deduced from InputVersion ("&
          & // i2c(implicitVersion) // ") differs from version explicitely set in ParserVersion ("&
          & // i2c(inputVersion) // ")", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      else
        inputVersion = implicitVersion
        call setChildValue(node, "ParserVersion", inputVersion, errStatus, child=child)
        @:PROPAGATE_ERROR(errStatus)
      end if
    else
      call getChildValue(node, "ParserVersion", inputVersion, errStatus, parserVersion, child=child)
      @:PROPAGATE_ERROR(errStatus)
    end if

    if (inputVersion < 1 .or. inputVersion > parserVersion) then
      call detailedError(child, "Invalid parser version (" // i2c(inputVersion) // ")", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else if (inputVersion < minVersion) then
      call detailedError(child, "Sorry, no compatibility mode for parser version "&
          & // i2c(inputVersion) // " (too old)", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else if (inputVersion /= parserVersion) then
      write(stdout, "(A,I2,A,I2,A)") "***  Converting input from parser version ",&
          & inputVersion, " to parser version ", parserVersion, " ..."
      call convertOldHSD(root, inputVersion, parserVersion, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      write(stdout, "(A,/)") "***  Done."
    end if

    call getChildValue(node, "WriteHSDInput", flags%tWriteHSD, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)
    if (.not. flags%tWriteHSD) then
      call detailedWarning(node, "WriteHSDInput turned off. You are not guaranteed" // newline // &
          &" to able to obtain the same results with a later version of the code!" // newline // &
          & "(the dftb_pin.hsd file DOES guarantee this)")
    end if
    call getChildValue(node, "StopAfterParsing", flags%tStop, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(node, "IgnoreUnprocessedNodes", flags%tIgnoreUnprocessed, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readParserOptions


  !> Read in Geometry
  subroutine readGeometry(node, input, errStatus)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Input structure to be filled
    type(TInputData), intent(inout) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child
    type(string) :: buffer

    call getChildValue(node, "", value1, errStatus, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    input%geom%tPeriodic = .false.
    input%geom%tHelical = .false.
    select case (char(buffer))
    case ("genformat")
      call readTGeometryGen(value1, input%geom, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("xyzformat")
      call readTGeometryXyz(value1, input%geom, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("vaspformat")
      call readTGeometryVasp(value1, input%geom, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("lammpsformat")
      call readTGeometryLammps(value1, input%geom, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case default
      call setUnprocessed(value1)
      call readTGeometryHSD(child, input%geom, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readGeometry


  !> Read in driver properties
#:if WITH_TRANSPORT
  subroutine readDriver(node, parent, geom, ctrl, transpar, errStatus)
#:else
  subroutine readDriver(node, parent, geom, ctrl, errStatus)
#:endif

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Parent of node (for error messages)
    type(fnode), pointer :: parent

    !> geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(in) :: transpar
  #:endif

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child, child2, child3, value1, value2, field

    type(string) :: buffer, buffer2, modifier
  #:if WITH_SOCKETS
    character(lc) :: sTmp
  #:endif

    ! Default range of atoms to move (may be adjusted if contacts present)
    character(mc) :: atomsRange

    character(mc) :: modeName
    logical :: isMaxStepNeeded

    ctrl%isGeoOpt = .false.
    ctrl%tCoordOpt = .false.
    ctrl%tLatOpt = .false.

    ctrl%iGeoOpt = geoOptTypes%none
    ctrl%tMD = .false.
    ctrl%iThermostat = 0
    ctrl%tForces = .false.
    ctrl%tSetFillingTemp = .false.

    atomsRange = "1:-1"
  #:if WITH_TRANSPORT
    if (transpar%defined) then
      ! only those atoms in the device region
      write(atomsRange,"(I0,':',I0)")transpar%idxdevice
    end if
  #:endif

    call renameChildren(parent, "GeometryOptimization", "GeometryOptimisation")
    call getNodeName2(node, buffer)
    driver: select case (char(buffer))
    case ("")
      modeName = ""
      continue
    case ("none")
      modeName = ""
      continue
    case ("geometryoptimisation")
      modeName = "geometry optimisation"

      if (geom%tHelical) then
        call detailedError(node, "GeometryOptimisation driver currently does not support helical&
            & geometries", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

      allocate(ctrl%geoOpt)

      call readGeoOptInput(node, geom, ctrl%geoOpt, atomsRange, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(node, "AppendGeometries", ctrl%tAppendGeo, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)

      ctrl%tForces = .true.
      ctrl%restartFreq = 1

    case ("steepestdescent")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ! Steepest downhill optimisation
      ctrl%iGeoOpt = geoOptTypes%steepestDesc
      call commonGeoOptions(node, ctrl, geom, atomsRange, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case ("conjugategradient")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')// "Please use the GeometryOptimisation driver instead.")

      ! Conjugate gradient location optimisation
      ctrl%iGeoOpt = geoOptTypes%conjugateGrad
      call commonGeoOptions(node, ctrl, geom, atomsRange, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case("gdiis")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ! Gradient DIIS optimisation, only stable in the quadratic region
      ctrl%iGeoOpt = geoOptTypes%diis
      call getChildValue(node, "alpha", ctrl%deltaGeoOpt, errStatus, 1.0E-1_dp)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(node, "Generations", ctrl%iGenGeoOpt, errStatus, 8)
      @:PROPAGATE_ERROR(errStatus)
      call commonGeoOptions(node, ctrl, geom, atomsRange, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case ("lbfgs")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ctrl%iGeoOpt = geoOptTypes%lbfgs

      allocate(ctrl%lbfgsInp)
      call getChildValue(node, "Memory", ctrl%lbfgsInp%memory, errStatus, 20)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(node, "LineSearch", ctrl%lbfgsInp%isLineSearch, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)

      isMaxStepNeeded = .not. ctrl%lbfgsInp%isLineSearch
      if (isMaxStepNeeded) then
        call getChildValue(node, "setMaxStep", ctrl%lbfgsInp%isLineSearch, errStatus,&
            & isMaxStepNeeded)
        @:PROPAGATE_ERROR(errStatus)
        ctrl%lbfgsInp%MaxQNStep = isMaxStepNeeded
      else
        call getChildValue(node, "oldLineSearch", ctrl%lbfgsInp%isOldLS, errStatus, .false.)
        @:PROPAGATE_ERROR(errStatus)
      end if

      call commonGeoOptions(node, ctrl, geom, atomsRange, errStatus,&
          & isMaxStepNeeded=ctrl%lbfgsInp%isLineSearch)
      @:PROPAGATE_ERROR(errStatus)

    case ("fire")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ctrl%iGeoOpt = geoOptTypes%fire
      call commonGeoOptions(node, ctrl, geom, atomsRange, errStatus, isMaxStepNeeded=.false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(node, "TimeStep", ctrl%deltaT, errStatus, 1.0_dp, modifier=modifier,&
          & child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), timeUnits, field, ctrl%deltaT, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case("secondderivatives")
      ! currently only numerical derivatives of forces is implemented

      modeName = "second derivatives"

      ctrl%tDerivs = .true.
      ctrl%tForces = .true.

      call getChildValue(node, "Atoms", buffer2, errStatus, trim(atomsRange), child=child,&
          & multiple=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getSelectedAtomIndices(child, char(buffer2), geom%speciesNames, geom%species,&
          & ctrl%indDerivAtom, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (size(ctrl%indDerivAtom) == 0) then
        @:RAISE_ERROR(errStatus, -1, "No atoms specified for derivatives calculation.")
      end if

      call getChild(node, "MovedAtoms", child, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child)) then
        if (.not. isContiguousRange(ctrl%indDerivAtom)) then
          call detailedError(child, "Atoms for calculation of partial Hessian must be a contiguous&
              & range.", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        call getChildValue(child, "", buffer2, errStatus, child=child2, multiple=.true.)
        @:PROPAGATE_ERROR(errStatus)
        call getSelectedAtomIndices(child2, char(buffer2), geom%speciesNames, geom%species,&
            & ctrl%indMovedAtom, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        if (.not. isContiguousRange(ctrl%indMovedAtom)) then
          call detailedError(child2, "MovedAtoms for calculation of partial Hessian must be a &
              & contiguous range.", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        if (.not. containsAll(ctrl%indDerivAtom, ctrl%indMovedAtom)) then
          call detailedError(child2, "MovedAtoms has indices not contained in Atoms.", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      else
        ctrl%indMovedAtom = ctrl%indDerivAtom
      end if
      ctrl%nrMoved = size(ctrl%indMovedAtom)

      call getChildValue(node, "Delta", ctrl%deriv2ndDelta, errStatus, 1.0E-4_dp,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, field, ctrl%deriv2ndDelta, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case ("velocityverlet")
      ! molecular dynamics

      modeName = "molecular dynamics"

      ctrl%tForces = .true.
      ctrl%tMD = .true.

      call getChildValue(node, "MDRestartFrequency", ctrl%restartFreq, errStatus, 1)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(node, "MovedAtoms", buffer2, errStatus, trim(atomsRange), child=child,&
          &multiple=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getSelectedAtomIndices(child, char(buffer2), geom%speciesNames, geom%species,&
          & ctrl%indMovedAtom, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      ctrl%nrMoved = size(ctrl%indMovedAtom)
      if (ctrl%nrMoved == 0) then
        @:RAISE_ERROR(errStatus, -1, "No atoms specified for molecular dynamics.")
      end if
      call readInitialVelocities(node, ctrl, geom%nAtom, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(node, "KeepStationary", ctrl%tMDstill, errStatus, .true.)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%tMDstill .and. geom%nAtom == 1) then
        @:RAISE_ERROR(errStatus, -1, "Removing translational freedom with only one atom not&
            & possible.")
      end if

      call getChildValue(node, "TimeStep", ctrl%deltaT, errStatus, modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), timeUnits, field, ctrl%deltaT, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(node, "Thermostat", value1, errStatus, child=child)
      @:PROPAGATE_ERROR(errStatus)
      call getNodeName(value1, buffer2)

      thermostat: select case(char(buffer2))
      case ("berendsen")
        ctrl%iThermostat = 2
        ! Read temperature or temperature profiles
        call getChildValue(value1, "Temperature", value2, errStatus, modifier=modifier,&
            & child=child2)
        @:PROPAGATE_ERROR(errStatus)
        call getNodeName(value2, buffer)

        select case(char(buffer))
        case (textNodeName)
          call readTemperature(child2, ctrl, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        case ("temperatureprofile")
          call readTemperatureProfile(value2, char(modifier), ctrl, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        case default
          call detailedError(value2, "Invalid method name.", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end select

        call getChild(value1, "CouplingStrength", child2, errStatus, requested=.false.)
        @:PROPAGATE_ERROR(errStatus)
        if (associated(child2)) then
          call getChildValue(child2, "", ctrl%wvScale, errStatus)
          @:PROPAGATE_ERROR(errStatus)
          call getChild(value1, "Timescale", child2, errStatus, modifier=modifier,&
              & requested=.false.)
          @:PROPAGATE_ERROR(errStatus)
          if (associated(child2)) then
            @:RAISE_ERROR(errStatus, -1, "Only Coupling strength OR Timescale can be set for&
                & Berendsen thermostats.")
          end if
        else
          call getChild(value1, "Timescale", child2, errStatus, modifier=modifier,&
              & requested=.false.)
          @:PROPAGATE_ERROR(errStatus)
          if (associated(child2)) then
            call getChildValue(child2, "", ctrl%wvScale, errStatus, modifier=modifier, child=child3)
            @:PROPAGATE_ERROR(errStatus)
            call convertUnitHsd(char(modifier), timeUnits, child3, ctrl%wvScale, errStatus)
            @:PROPAGATE_ERROR(errStatus)
            ctrl%wvScale = ctrl%deltaT / ctrl%wvScale
          else
            @:RAISE_ERROR(errStatus, -1, "Either CouplingStrength or Timescale must be set&
                & for Berendsen thermostats.")
          end if
        end if

      case ("nosehoover")
        ctrl%iThermostat = 3
        ! Read temperature or temperature profiles
        call getChildValue(value1, "Temperature", value2, errStatus, modifier=modifier,&
            & child=child2)
        @:PROPAGATE_ERROR(errStatus)
        call getNodeName(value2, buffer)

        select case(char(buffer))
        case (textNodeName)
          call readTemperature(child2, ctrl, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        case ("temperatureprofile")
          call readTemperatureProfile(value2, char(modifier), ctrl, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        case default
          call detailedError(value2, "Invalid method name.", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end select

        call getChildValue(value1, "CouplingStrength", ctrl%wvScale, errStatus, modifier=modifier,&
            & child=field)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(modifier), freqUnits, field, ctrl%wvScale, errStatus)
        @:PROPAGATE_ERROR(errStatus)

        call getChildValue(value1, "ChainLength", ctrl%nh_npart, errStatus, 3)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(value1, "Order", ctrl%nh_nys, errStatus, 3)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(value1, "IntegratorSteps", ctrl%nh_nc, errStatus, 1)
        @:PROPAGATE_ERROR(errStatus)

        call getChild(value1, "Restart", child3, errStatus, requested=.false.)
        @:PROPAGATE_ERROR(errStatus)
        if (associated(child3)) then
          allocate(ctrl%xnose(ctrl%nh_npart))
          allocate(ctrl%vnose(ctrl%nh_npart))
          allocate(ctrl%gnose(ctrl%nh_npart))
          call getChildValue(child3, "x", ctrl%xnose, errStatus)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(child3, "v", ctrl%vnose, errStatus)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(child3, "g", ctrl%gnose, errStatus)
          @:PROPAGATE_ERROR(errStatus)
          ctrl%tInitNHC = .true.
        else
          ctrl%tInitNHC = .false.
        end if

      case ("andersen")
        ctrl%iThermostat = 1
        ! Read temperature or temperature profiles
        call getChildValue(value1, "Temperature", value2, errStatus, modifier=modifier,&
            & child=child2)
        @:PROPAGATE_ERROR(errStatus)
        call getNodeName(value2, buffer)

        select case(char(buffer))
        case (textNodeName)
          call readTemperature(child2, ctrl, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        case ("temperatureprofile")
          call readTemperatureProfile(value2, char(modifier), ctrl, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        case default
          call detailedError(value2, "Invalid method name.", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end select

        call getChildValue(value1, "ReselectProbability", ctrl%wvScale, errStatus, child=child3)
        @:PROPAGATE_ERROR(errStatus)
        if (ctrl%wvScale <= 0.0_dp .or. ctrl%wvScale > 1.0_dp) then
          call detailedError(child3, "ReselectProbability must be in the range (0,1]!", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        call getChildValue(value1, "ReselectIndividually", ctrl%tRescale, errStatus)
        @:PROPAGATE_ERROR(errStatus)

      case ("none")
        ctrl%iThermostat = 0
        allocate(ctrl%tempSteps(1))
        allocate(ctrl%tempValues(1))

        if (ctrl%tReadMDVelocities) then
          ! without a thermostat, if we know the initial velocities, we do not
          ! need a temperature, so just set it to something 'safe'
          ctrl%tempAtom = minTemp
        else
          call readMDInitTemp(value1, ctrl%tempAtom, minTemp, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      case default
        call getNodeHSDName(value1, buffer2)
        call detailedError(child, "Invalid thermostat '" // char(buffer2) // "'", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end select thermostat

      if (ctrl%maxRun < -1) then
        call getChildValue(node, "Steps", ctrl%maxRun, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

      call getChildValue(node, "OutputPrefix", buffer2, errStatus, "geo_end")
      @:PROPAGATE_ERROR(errStatus)
      ctrl%outFile = unquote(char(buffer2))

      call getChildValue(node, "Plumed", ctrl%tPlumed, errStatus, default=.false., child=child)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%tPlumed .and. .not. withPlumed) then
        call detailedError(child, "Metadynamics can not be used since code has been compiled&
            & without PLUMED support", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

      if (geom%tPeriodic) then

        call getChild(node, "Barostat", child, errStatus, requested=.false.)
        @:PROPAGATE_ERROR(errStatus)
        if (.not. associated(child)) then
          call setChild(node, "Barostat", child, errStatus)
          @:PROPAGATE_ERROR(errStatus)
          ctrl%tBarostat = .false.
          ctrl%pressure = 0.0_dp
          ctrl%BarostatStrength = 0.0_dp
        else
          if (ctrl%nrMoved /= geom%nAtom) then
            @:RAISE_ERROR(errStatus, -1, "Dynamics for a subset of atoms is not currently&
                & possible when using a barostat")
          end if
          call getChildValue(child, "Pressure", ctrl%pressure, errStatus, modifier=modifier,&
              & child=child2)
          @:PROPAGATE_ERROR(errStatus)
          call convertUnitHsd(char(modifier), pressureUnits, child2, ctrl%pressure, errStatus)
          @:PROPAGATE_ERROR(errStatus)
          call getChild(child, "Coupling", child2, errStatus, requested=.false.)
          @:PROPAGATE_ERROR(errStatus)
          if (associated(child2)) then
            call getChildValue(child2, "", ctrl%BarostatStrength, errStatus)
            @:PROPAGATE_ERROR(errStatus)
            call getChild(child, "Timescale", child2, errStatus, modifier=modifier,&
                & requested=.false.)
            @:PROPAGATE_ERROR(errStatus)
            if (associated(child2)) then
              @:RAISE_ERROR(errStatus, -1, "Only Coupling strength OR Timescale can be set for&
                  & Barostatting.")
            end if
          else
            call getChild(child, "Timescale", child2, errStatus, modifier=modifier,&
                & requested=.false.)
            @:PROPAGATE_ERROR(errStatus)
            if (associated(child2)) then
              call getChildValue(child2, "", ctrl%BarostatStrength, errStatus, modifier=modifier,&
                  & child=child3)
              @:PROPAGATE_ERROR(errStatus)
              call convertUnitHsd(char(modifier), timeUnits, child3, ctrl%BarostatStrength,&
                  & errStatus)
              @:PROPAGATE_ERROR(errStatus)
              ctrl%BarostatStrength = ctrl%deltaT / ctrl%BarostatStrength
            else
              @:RAISE_ERROR(errStatus, -1, "Either Coupling strength or Timescale must be set&
                  & for Barostatting.")
            end if
          end if
          call getChildValue(child, "Isotropic", ctrl%tIsotropic, errStatus, .true.)
          @:PROPAGATE_ERROR(errStatus)
          ctrl%tBarostat = .true.
        end if
      end if

      if (ctrl%hamiltonian == hamiltonianTypes%dftb) then
        call readXlbomdOptions(node, ctrl%xlbomd, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

      call getInputMasses(node, geom, ctrl%masses, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case ("socket")
      ! external socket control of the run (once initialised from input)

      modeName = "socket control"

    #:if WITH_SOCKETS
      ctrl%tForces = .true.
      allocate(ctrl%socketInput)
      call getChild(node, 'File', child2, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChild(node, 'Host', child3, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child2) .eqv. associated(child3)) then
        @:RAISE_ERROR(errStatus, -1, 'Either Host or File (but not both) must be set for socket&
            & communication')
      end if

      ! File communication
      if (associated(child2)) then
        call getChildValue(child2, "", buffer2, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        ctrl%socketInput%host = unquote(char(buffer2))
        ! zero it to signal to initprogram to use a unix file
        ctrl%socketInput%port = 0
      else
        call getChildValue(child3, "", buffer2, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        ctrl%socketInput%host = unquote(char(buffer2))
        call getChildValue(node, "Port", ctrl%socketInput%port, errStatus, child=field)
        @:PROPAGATE_ERROR(errStatus)
        if (ctrl%socketInput%port <= 0) then
          call detailedError(field, "Invalid port number", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      end if

      call getChildValue(node, "Protocol", value1, errStatus, "i-PI", child=child)
      @:PROPAGATE_ERROR(errStatus)
      call getNodeName(value1, buffer)
      select case(char(buffer))
      case("i-pi")
        ctrl%socketInput%protocol = IPI_PROTOCOLS%IPI_1
        ! want a file path
        if (ctrl%socketInput%port == 0) then
          call getChildValue(node, "Prefix", buffer2, errStatus, "/tmp/ipi_")
          @:PROPAGATE_ERROR(errStatus)
          sTmp = unquote(char(buffer2))
          ctrl%socketInput%host = trim(sTmp) // trim(ctrl%socketInput%host)
        end if

      case default
        call detailedError(child, "Invalid protocol '" // char(buffer) // "'", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end select
      call getChildValue(node, "Verbosity", ctrl%socketInput%verbosity, errStatus, 0)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(node, "MaxSteps", ctrl%maxRun, errStatus, 200)
      @:PROPAGATE_ERROR(errStatus)

    #:else
      call detailedError(node, "Program had been compiled without socket support", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:endif

    case default

      call getNodeHSDName(node, buffer)
      call detailedError(parent, "Invalid driver '" // char(buffer) // "'", errStatus)
      @:PROPAGATE_ERROR(errStatus)

    end select driver

  #:if WITH_TRANSPORT
    if (ctrl%solver%isolver == electronicSolverTypes%OnlyTransport .and. trim(modeName) /= "") then
      call detailederror(node, "transportOnly solver cannot be used with "//trim(modeName),&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
  #:endif

  end subroutine readDriver


  !> Simple function to check that an array of indices is a contigous range
  function isContiguousRange(indices) result(isContiguous)

    !> Array of atomic indices
    integer, intent(in) :: indices(:)

    !> whether indices are contigous
    logical :: isContiguous

    isContiguous = all(indices(: size(indices) - 1) + 1 == indices(2:))

  end function isContiguousRange


  !> checks that the array subindices is contained in indices
  function containsAll(indices, subindices)

    !> Array of atomic indices to check against
    integer, intent(in) :: indices(:)

    !> Array of atomic indices to check
    integer, intent(in) :: subindices(:)

    !> whether indices are contigous
    logical :: containsAll

    integer :: kk

    containsAll = .false.
    do kk = 1, size(subindices)
      if (.not. any(indices == subindices(kk))) return
    end do
    containsAll = .true.

  end function containsAll


  !> Common geometry optimisation settings for various drivers
  subroutine commonGeoOptions(node, ctrl, geom, atomsRange, errStatus, isMaxStepNeeded)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Default range of moving atoms (may be restricted for example by contacts in transport
    !> calculations)
    character(len=*), intent(in) :: atomsRange

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> Is the maximum step size relevant for this driver
    logical, intent(in), optional :: isMaxStepNeeded

    type(fnode), pointer :: child, field
    type(string) :: buffer2, modifier
    logical :: isMaxStep

    if (present(isMaxStepNeeded)) then
      isMaxStep = isMaxStepNeeded
    else
      isMaxStep = .true.
    end if

    ctrl%tForces = .true.
    ctrl%restartFreq = 1

    call getChildValue(node, "LatticeOpt", ctrl%tLatOpt, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    if (ctrl%tLatOpt) then
      call getChildValue(node, "Pressure", ctrl%pressure, errStatus, 0.0_dp, modifier=modifier,&
          & child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), pressureUnits, child, ctrl%pressure, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(node, "FixAngles", ctrl%tLatOptFixAng, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%tLatOptFixAng) then
        call getChildValue(node, "FixLengths", ctrl%tLatOptFixLen, errStatus,&
            & [.false.,.false.,.false.])
        @:PROPAGATE_ERROR(errStatus)
      else
        call getChildValue(node, "Isotropic", ctrl%tLatOptIsotropic, errStatus, .false.)
        @:PROPAGATE_ERROR(errStatus)
      end if
      if (isMaxStep) then
        call getChildValue(node, "MaxLatticeStep", ctrl%maxLatDisp, errStatus, 0.2_dp)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end if
    call getChildValue(node, "MovedAtoms", buffer2, errStatus, trim(atomsRange), child=child,&
        & multiple=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call getSelectedAtomIndices(child, char(buffer2), geom%speciesNames, geom%species,&
        & ctrl%indMovedAtom, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ctrl%nrMoved = size(ctrl%indMovedAtom)
    ctrl%tCoordOpt = (ctrl%nrMoved /= 0)
    if (ctrl%tCoordOpt) then
      if (isMaxStep) then
        call getChildValue(node, "MaxAtomStep", ctrl%maxAtomDisp, errStatus, 0.2_dp)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end if
    call getChildValue(node, "MaxForceComponent", ctrl%maxForce, errStatus, 1e-4_dp,&
        & modifier=modifier, child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), forceUnits, field, ctrl%maxForce, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "MaxSteps", ctrl%maxRun, errStatus, 200)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "StepSize", ctrl%deltaT, errStatus, 100.0_dp, modifier=modifier,&
        & child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), timeUnits, field, ctrl%deltaT, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "OutputPrefix", buffer2, errStatus, "geo_end")
    @:PROPAGATE_ERROR(errStatus)
    ctrl%outFile = unquote(char(buffer2))
    call getChildValue(node, "AppendGeometries", ctrl%tAppendGeo, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call readGeoConstraints(node, ctrl, geom%nAtom, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (ctrl%tLatOpt) then
      if (ctrl%nrConstr/=0) then
        @:RAISE_ERROR(errStatus, -1, "Lattice optimisation and constraints currently incompatible.")
      end if
      if (ctrl%nrMoved/=0.and.ctrl%nrMoved<geom%nAtom) then
        @:RAISE_ERROR(errStatus, -1, "Subset of optimising atoms not currently possible with&
            & lattice optimisation.")
      end if
    end if
    ctrl%isGeoOpt = ctrl%tLatOpt .or. ctrl%tCoordOpt

  end subroutine commonGeoOptions


  !> Extended lagrangian options for XLBOMD
  subroutine readXlbomdOptions(node, input, errStatus)

    !> node in the input tree
    type(fnode), pointer :: node

    !> extracted settings on exit
    type(TXLBOMDInp), allocatable, intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: pXlbomd, pXlbomdFast, pRoot, pChild
    integer :: nKappa
    logical :: tXlbomdFast

    call getChild(node, 'Xlbomd', pXlbomd, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChild(node, 'XlbomdFast', pXlbomdFast, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (.not. (associated(pXlbomd) .or. associated(pXlbomdFast))) then
      return
    end if
    if (associated(pXlbomd) .and. associated(pXlbomdFast)) then
      call detailedError(pXlbomdFast, "Blocks 'Xlbomd' and 'XlbomdFast' are mutually exclusive",&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    if (associated(pXlbomdFast)) then
      tXlbomdFast = .true.
      pRoot => pXlbomdFast
    else
      tXlbomdFast = .false.
      pRoot => pXlbomd
    end if
    allocate(input)
    call getChildValue(pRoot, 'IntegrationSteps', input%nKappa, errStatus, 5, child=pChild)
    @:PROPAGATE_ERROR(errStatus)
    ! Workaround for nagfor 6.1 as the derived type in the array comparison
    ! is not recognized as such
    nKappa = input%nKappa
    if (all([5, 6, 7] /= nKappa)) then
      call detailedError(pChild, 'Invalid number of integration steps (must be 5, 6 or 7)',&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    call getChildValue(pRoot, 'PreSteps', input%nPreSteps, errStatus, 0)
    @:PROPAGATE_ERROR(errStatus)

    ! Since inverse Jacobian is not enabled, we can set FullSccSteps
    ! to its minimal value (no averaging of inverse Jacobians is done)
    !call getChildValue(child, 'FullSccSteps', input%nFullSccSteps, errStatus,&
    !    & input%nKappa + 1, child=child2)
    ! @:PROPAGATE_ERROR(errStatus)
    input%nFullSccSteps = input%nKappa + 1
    !if (input%nFullSccSteps < input%nKappa + 1) then
    !  call detailedError(child2, 'Nr. of full SCC steps must be greater by&
    !      & one than integration steps', errStatus)
    ! @:PROPAGATE_ERROR(errStatus)
    !end if

    if (tXlbomdFast) then
      call getChildValue(pRoot, 'TransientSteps', input%nTransientSteps, errStatus, 10)
      @:PROPAGATE_ERROR(errStatus)
      input%minSccIter = 1
      input%maxSccIter = 1
      ! Dummy value as minSccIter and maxSccIter have been set to 1.
      input%sccTol = 1e-5_dp
      call getChildValue(pRoot, 'Scale', input%scale, errStatus, 1.0_dp, child=pChild)
      @:PROPAGATE_ERROR(errStatus)
      if (input%scale <= 0.0_dp .or. input%scale > 1.0_dp) then
        call detailedError(pChild, 'Scaling value must be in the interval (0.0, 1.0]', errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

      !! Inverse Jacobian is experimental feature so far.
      !call getChildValue(child, "UseJacobianKernel", &
      !    & input%useInverseJacobian, errStatus, .false.)
      ! @:PROPAGATE_ERROR(errStatus)
      !if (input%useInverseJacobian) then
      !  call getChildValue(child, "ReadJacobianKernel", &
      !      & input%readInverseJacobian, errStatus, .false.)
      !  @:PROPAGATE_ERROR(errStatus)
      !end if
      input%useInverseJacobian = .false.
      input%readInverseJacobian = .false.
    else
      input%nTransientSteps = 0
      call getChildValue(pRoot, 'MinSccIterations', input%minSCCIter, errStatus, 1)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(pRoot, 'MaxSccIterations', input%maxSCCIter, errStatus, 200)
      @:PROPAGATE_ERROR(errStatus)
      if (input%maxSCCIter <= 0) then
        call detailedError(pRoot,"MaxSccIterations must be >= 1", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getChildValue(pRoot, 'SccTolerance', input%sccTol, errStatus, 1e-5_dp)
      @:PROPAGATE_ERROR(errStatus)
      input%scale = 1.0_dp
      input%useInverseJacobian = .false.
      input%readInverseJacobian = .false.
    end if

  end subroutine readXlbomdOptions


  !> Reads geometry constraints
  subroutine readGeoConstraints(node, ctrl, nAtom, errStatus)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child
    type(string) :: buffer
    type(TListIntR1) :: intBuffer
    type(TListRealR1) :: realBuffer

    call getChildValue(node, "Constraints", value1, errStatus, "", child=child,&
        & allowEmptyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName2(value1, buffer)
    if (char(buffer) == "") then
      ctrl%nrConstr = 0
    else
      call init(intBuffer)
      call init(realBuffer)
      call getChildValue(child, "", 1, intBuffer, 3, realBuffer, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      ctrl%nrConstr = len(intBuffer)
      allocate(ctrl%conAtom(ctrl%nrConstr))
      allocate(ctrl%conVec(3, ctrl%nrConstr))
      call asVector(intBuffer, ctrl%conAtom)
      if (.not.all(ctrl%conAtom<=nAtom)) then
        call detailedError(node, "Non-existent atom specified in constraint", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call asArray(realBuffer, ctrl%conVec)
      call destruct(intBuffer)
      call destruct(realBuffer)
    end if

  end subroutine readGeoConstraints


  !> Reads MD velocities
  subroutine readInitialVelocities(node, ctrl, nAtom, errStatus)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Total number of all atoms
    integer, intent(in) :: nAtom

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child
    type(string) :: buffer, modifier
    type(TListRealR1) :: realBuffer
    integer :: nVelocities
    real(dp), allocatable :: tmpVelocities(:,:)

    call getChildValue(node, "Velocities", value1, errStatus, "", child=child, modifier=modifier,&
        & allowEmptyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName2(value1, buffer)
    if (char(buffer) == "") then
      ctrl%tReadMDVelocities = .false.
    else
      call init(realBuffer)
      call getChildValue(child, "", 3, realBuffer, errStatus, modifier=modifier)
      @:PROPAGATE_ERROR(errStatus)
      nVelocities = len(realBuffer)
      if (nVelocities /= nAtom) then
        call detailedError(node, "Incorrect number of specified velocities: " // i2c(3*nVelocities)&
            & // " supplied, " // i2c(3*nAtom) // " required.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      allocate(tmpVelocities(3, nVelocities))
      call asArray(realBuffer, tmpVelocities)
      if (len(modifier) > 0) then
        call convertUnitHsd(char(modifier), VelocityUnits, child, tmpVelocities, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call destruct(realBuffer)
      allocate(ctrl%initialVelocities(3, ctrl%nrMoved))
      ctrl%initialVelocities(:,:) = tmpVelocities(:,ctrl%indMovedAtom(:))
      ctrl%tReadMDVelocities = .true.
    end if

  end subroutine readInitialVelocities


  !> Reads atomic masses from input file, eventually overwriting those in the SK files
  subroutine getInputMasses(node, geo, masses, errStatus)

    !> relevant node of input data
    type(fnode), pointer :: node

    !> geometry object, which contains atomic species information
    type(TGeometry), intent(in) :: geo

    !> masses to be returned
    real(dp), allocatable, intent(out) :: masses(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child, child2, child3, val
    type(fnodeList), pointer :: children
    integer, allocatable :: pTmpI1(:)
    type(string) :: buffer, modifier
    real(dp) :: rTmp
    integer :: ii, jj, iAt

    call getChildValue(node, "Masses", val, errStatus, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)
    @:PROPAGATE_ERROR(errStatus)

    ! Read individual atom specifications
    call getChildren(child, "Mass", children)
    if (getLength(children) == 0) then
      call destroyNodeList(children)
      return
    end if

    allocate(masses(geo%nAtom))
    masses(:) = -1.0_dp
    do ii = 1, getLength(children)
      call getItem1(children, ii, child2)
      call getChildValue(child2, "Atoms", buffer, errStatus, child=child3, multiple=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species, pTmpI1,&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child2, "MassPerAtom", rTmp, errStatus, modifier=modifier, child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), massUnits, child, rTmp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      do jj = 1, size(pTmpI1)
        iAt = pTmpI1(jj)
        if (masses(iAt) >= 0.0_dp) then
          call detailedWarning(child3, "Previous setting for the mass  of atom" // i2c(iAt) //&
              & " overwritten")
        end if
        masses(iAt) = rTmp
      end do
      deallocate(pTmpI1)
    end do
    call destroyNodeList(children)

  end subroutine getInputMasses


  !> Reads Hamiltonian
#:if WITH_TRANSPORT
  subroutine readHamiltonian(node, ctrl, geo, slako, tp, greendens, poisson, errStatus)
#:else
  subroutine readHamiltonian(node, ctrl, geo, slako, poisson, errStatus)
#:endif

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Slater-Koster structure to be filled
    type(TSlater), intent(inout) :: slako

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens
  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(string) :: buffer

    call getNodeName(node, buffer)
    select case (char(buffer))
    case ("dftb")
  #:if WITH_TRANSPORT
      call readDFTBHam(node, ctrl, geo, slako, tp, greendens, poisson, errStatus)
      @:PROPAGATE_ERROR(errStatus)
  #:else
      call readDFTBHam(node, ctrl, geo, slako, poisson, errStatus)
      @:PROPAGATE_ERROR(errStatus)
  #:endif
    case ("xtb")
  #:if WITH_TRANSPORT
      call readXTBHam(node, ctrl, geo, tp, greendens, poisson, errStatus)
      @:PROPAGATE_ERROR(errStatus)
  #:else
      call readXTBHam(node, ctrl, geo, poisson, errStatus)
      @:PROPAGATE_ERROR(errStatus)
  #:endif
    case default
      call detailedError(node, "Invalid Hamiltonian", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readHamiltonian


  !> Reads DFTB-Hamiltonian
#:if WITH_TRANSPORT
  subroutine readDFTBHam(node, ctrl, geo, slako, tp, greendens, poisson, errStatus)
#:else
  subroutine readDFTBHam(node, ctrl, geo, slako, poisson, errStatus)
#:endif

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Slater-Koster structure to be filled
    type(TSlater), intent(inout) :: slako

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child, child2, child3
    type(fnodeList), pointer :: children
    type(string) :: buffer, buffer2, modifier
    type(TListInt) :: li
    type(TListInt), allocatable :: liN(:)
    type(TListIntR1), allocatable :: li1N(:)
    type(TListReal), allocatable :: lrN(:)
    type(TListCharLc), allocatable :: skFiles(:,:)
    type(TListString) :: lStr
    type(TListIntR1), allocatable :: angShells(:)
    logical, allocatable :: repPoly(:,:)
    integer :: iSp1, iSp2, ii
    character(lc) :: prefix, suffix, separator, elem1, elem2, strTmp
    character(lc) :: errorStr
    logical :: tLower, tExist
    integer, allocatable :: pTmpI1(:)
    real(dp) :: rTmp
    integer, allocatable :: iTmpN(:)
    integer :: nShell, skInterMeth
    real(dp) :: rSKCutOff
    type(string), allocatable :: searchPath(:)
    character(len=:), allocatable :: strOut

    !> For range separation
    type(TRangeSepSKTag) :: rangeSepSK

    ctrl%hamiltonian = hamiltonianTypes%dftb

    call readMaxAngularMomentum(node, geo, angShells, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Orbitals and angular momenta for the given shells (once the SK files contain the full
    ! information about the basis, this will be moved to the SK reading routine).
    allocate(slako%orb)
    call setupOrbitals(slako%orb, geo, angShells, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Slater-Koster files
    call getParamSearchPath(searchPath)
    allocate(skFiles(geo%nSpecies, geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      do iSp2 = 1, geo%nSpecies
        call init(skFiles(iSp2, iSp1))
      end do
    end do
    call getChildValue(node, "SlaterKosterFiles", value1, errStatus, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case ("type2filenames")
      call getChildValue(value1, "Prefix", buffer2, errStatus, "")
      @:PROPAGATE_ERROR(errStatus)
      prefix = unquote(char(buffer2))
      call getChildValue(value1, "Suffix", buffer2, errStatus, "")
      @:PROPAGATE_ERROR(errStatus)
      suffix = unquote(char(buffer2))
      call getChildValue(value1, "Separator", buffer2, errStatus, "")
      @:PROPAGATE_ERROR(errStatus)
      separator = unquote(char(buffer2))
      call getChildValue(value1, "LowerCaseTypeName", tLower, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
      do iSp1 = 1, geo%nSpecies
        if (tLower) then
          elem1 = tolower(geo%speciesNames(iSp1))
        else
          elem1 = geo%speciesNames(iSp1)
        end if
        do iSp2 = 1, geo%nSpecies
          if (tLower) then
            elem2 = tolower(geo%speciesNames(iSp2))
          else
            elem2 = geo%speciesNames(iSp2)
          end if
          strTmp = trim(prefix) // trim(elem1) // trim(separator) &
              &// trim(elem2) // trim(suffix)
          call findFile(searchPath, strTmp, strOut)
          if (allocated(strOut)) strTmp = strOut
          call append(skFiles(iSp2, iSp1), strTmp)
          inquire(file=strTmp, exist=tExist)
          if (.not. tExist) then
            call detailedError(value1, "SK file with generated name '" // trim(strTmp)&
                & // "' does not exist.", errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end if
        end do
      end do
    case default
      call setUnprocessed(value1)
      do iSp1 = 1, geo%nSpecies
        do iSp2 = 1, geo%nSpecies
          strTmp = trim(geo%speciesNames(iSp1)) // "-" &
              &// trim(geo%speciesNames(iSp2))
          call init(lStr)
          call getChildValue(child, trim(strTmp), lStr, errStatus, child=child2)
          @:PROPAGATE_ERROR(errStatus)
          if (len(lStr) /= len(angShells(iSp1)) * len(angShells(iSp2))) then
            call detailedError(child2, "Incorrect number of Slater-Koster files", errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end if
          do ii = 1, len(lStr)
            call get(lStr, strTmp, ii)
            call findFile(searchPath, strTmp, strOut)
            if (allocated(strOut)) strTmp = strOut
            inquire(file=strTmp, exist=tExist)
            if (.not. tExist) then
              call detailedError(child2, "SK file '" // trim(strTmp) // "' does not exist'",&
                  & errStatus)
              @:PROPAGATE_ERROR(errStatus)
            end if
            call append(skFiles(iSp2, iSp1), strTmp)
          end do
          call destruct(lStr)
        end do
      end do
    end select

    ! Which repulsive is defined by polynomial? (Default: None)
    allocate(repPoly(geo%nSpecies, geo%nSpecies))
    call getChildValue(node, "PolynomialRepulsive", value1, errStatus, "", child=child,&
        & list=.true., allowEmptyValue=.true., dummyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName2(value1, buffer)
    select case (char(buffer))
    case ("")
      repPoly(:,:) = .false.
    case("setforall")
      call getChildValue(value1, "", repPoly(1,1), errStatus)
      @:PROPAGATE_ERROR(errStatus)
      repPoly(:,:) = repPoly(1,1)
    case default
      do iSp1 = 1, geo%nSpecies
        do iSp2 = 1, geo%nSpecies
          strTmp = trim(geo%speciesNames(iSp1)) // "-" &
              &// trim(geo%speciesNames(iSp2))
          call getChildValue(child, trim(strTmp), repPoly(iSp2, iSp1), errStatus, .false.)
          @:PROPAGATE_ERROR(errStatus)
        end do
      end do
      if (.not. all(repPoly .eqv. transpose(repPoly))) then
        call detailedError(value1, "Asymmetric definition (both A-B and B-A must&
            & be defined for using polynomial repulsive)", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end select

    call parseChimes(node, ctrl%chimesRepInput, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! SCC
    call getChildValue(node, "SCC", ctrl%tSCC, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)

    if (ctrl%tSCC) then
      call getChildValue(node, "ShellResolvedSCC", ctrl%tShellResolved, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
    else
      ctrl%tShellResolved = .false.
    end if

    call getChildValue(node, "OldSKInterpolation", ctrl%oldSKInter, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    if (ctrl%oldSKInter) then
      skInterMeth = skEqGridOld
    else
      skInterMeth = skEqGridNew
    end if

    call parseRangeSeparated(node, ctrl%rangeSepInp, skFiles, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if (.not. allocated(ctrl%rangeSepInp)) then
      call getChild(node, "TruncateSKRange", child, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child)) then
        call warning("Artificially truncating the SK table, this is normally a bad idea!")
        call SKTruncations(child, rSKCutOff, skInterMeth, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
            & skInterMeth, repPoly, errStatus, rSKCutOff)
        @:PROPAGATE_ERROR(errStatus)
      else
        rSKCutOff = 0.0_dp
        call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
            & skInterMeth, repPoly, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

    else

      call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
          & skInterMeth, repPoly, errStatus, rangeSepSK=rangeSepSK)
      @:PROPAGATE_ERROR(errStatus)
      ctrl%rangeSepInp%omega = rangeSepSk%omega
    end if


    do iSp1 = 1, geo%nSpecies
      call destruct(angShells(iSp1))
      do iSp2 = 1, geo%nSpecies
        call destruct(skFiles(iSp2, iSp1))
      end do
    end do
    deallocate(angShells)
    deallocate(skFiles)
    deallocate(repPoly)

    ! SCC parameters
    ifSCC: if (ctrl%tSCC) then

      ! get charge mixing options
      call readSccOptions(node, ctrl, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      ! DFTB hydrogen bond corrections
      call readHCorrection(node, geo, ctrl, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      !> TI-DFTB varibles for Delta DFTB
      call getChild(node, "NonAufbau", child, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child)) then
        ctrl%isNonAufbau = .true.
        call getChildValue(child, "SpinPurify", ctrl%isSpinPurify, errStatus, .true.)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child, "GroundGuess", ctrl%isGroundGuess, errStatus, .false.)
        @:PROPAGATE_ERROR(errStatus)
        ctrl%tSpin = .true.
        ctrl%t2Component = .false.
        ctrl%nrSpinPol = 0.0_dp
        ctrl%tSpinSharedEf = .false.
      else
        ctrl%isNonAufbau = .false.
      end if

    end if ifSCC

    ! Customize the reference atomic charges for virtual doping
    call readCustomReferenceOcc(node, slako%orb, slako%skOcc, geo, ctrl%customOccAtoms,&
        & ctrl%customOccFillings, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Spin calculation
    if (ctrl%reksInp%reksAlg == reksTypes%noReks  .and. .not.ctrl%isNonAufbau) then
    #:if WITH_TRANSPORT
      call readSpinPolarisation(node, ctrl, geo, tp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:else
      call readSpinPolarisation(node, ctrl, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:endif
    end if

    ! temporararily removed until debugged
    !if (.not. ctrl%tscc) then
    !  !! In a non-SCC calculation it is possible to upload charge shifts
    !  !! This is useful if the calculation can jump directly to the Analysis block
    !  call getChildValue(node, "ReadShifts", ctrl%tReadShifts, errStatus, .false.)
    ! @:PROPAGATE_ERROR(errStatus)
    !end if
    ctrl%tReadShifts = .false.

    ! External fields and potentials
    call readExternal(node, ctrl, geo, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Non-self-consistent spin-orbit coupling
    call readSpinOrbit(node, ctrl, geo, slako%orb, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Electronic solver
  #:if WITH_TRANSPORT
    call readSolver(node, ctrl, geo, tp, greendens, poisson, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if (tp%taskUpload) then
      ! Initialise variable, but unused
      ctrl%nrChrg =  0.0_dp
    else
      ! Charge
      call getChildValue(node, "Charge", ctrl%nrChrg, errStatus, 0.0_dp)
      @:PROPAGATE_ERROR(errStatus)
    end if
  #:else
    call readSolver(node, ctrl, geo, poisson, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Charge
    call getChildValue(node, "Charge", ctrl%nrChrg, errStatus, 0.0_dp)
    @:PROPAGATE_ERROR(errStatus)
  #:endif

    ! K-Points
    call readKPoints(node, ctrl, geo, ctrl%poorKSampling, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChild(node, "OrbitalPotential", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      allocate(ctrl%dftbUInp)
      call getChildValue(child, "Functional", buffer, errStatus, "fll")
      @:PROPAGATE_ERROR(errStatus)
      select case(tolower(char(buffer)))
      case ("fll")
        ctrl%dftbUInp%iFunctional = plusUFunctionals%fll
      case ("psic")
        ctrl%dftbUInp%iFunctional = plusUFunctionals%pSic
      case default
        call detailedError(child,"Unknown orbital functional :"// char(buffer), errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end select

      allocate(ctrl%dftbUInp%nUJ(geo%nSpecies))
      ctrl%dftbUInp%nUJ(:) = 0

      ! to hold list of U-J values for each atom
      allocate(lrN(geo%nSpecies))
      ! to hold count of U-J values for each atom
      allocate(liN(geo%nSpecies))
      ! to hold list of shells for each U-J block of values
      allocate(li1N(geo%nSpecies))

      do iSp1 = 1, geo%nSpecies
        call init(lrN(iSp1))
        call init(liN(iSp1))
        call init(li1N(iSp1))
        call getChildren(child, trim(geo%speciesNames(iSp1)), children)
        ctrl%dftbUInp%nUJ(iSp1) = getLength(children)
        do ii = 1, ctrl%dftbUInp%nUJ(iSp1)
          call getItem1(children, ii, child2)

          call init(li)
          call getChildValue(child2, "Shells", li, errStatus)
          @:PROPAGATE_ERROR(errStatus)
          allocate(pTmpI1(len(li)))
          call asArray(li,pTmpI1)
          call append(li1N(iSp1),pTmpI1)
          call append(liN(iSp1),size(pTmpI1))
          deallocate(pTmpI1)
          call destruct(li)
          call getChildValue(child2, "uj", rTmp, errStatus, 0.0_dp, modifier=modifier, child=child3)
          @:PROPAGATE_ERROR(errStatus)
          call convertUnitHsd(char(modifier), energyUnits, child3, rTmp, errStatus)
          @:PROPAGATE_ERROR(errStatus)
          if (rTmp < 0.0_dp) then
            write(errorStr,"(F12.8)")rTmp
            call detailedError(child2,"Negative value of U-J:" // errorStr, errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end if
          if (rTmp <= 1.0E-10_dp) then
            write(errorStr,"(F12.8)")rTmp
            call detailedError(child2,"Invalid value of U-J, too small: " // errorStr, errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end if
          call append(lrN(iSp1),rTmp)
        end do
        call destroyNodeList(children)
      end do

      do iSp1 = 1, geo%nSpecies
        ctrl%dftbUInp%nUJ(iSp1) = len(lrN(iSp1))
      end do
      allocate(ctrl%dftbUInp%UJ(maxval(ctrl%dftbUInp%nUJ),geo%nSpecies))
      ctrl%dftbUInp%UJ(:,:) = 0.0_dp
      allocate(ctrl%dftbUInp%niUJ(maxval(ctrl%dftbUInp%nUJ),geo%nSpecies))
      ctrl%dftbUInp%niUJ(:,:) = 0
      do iSp1 = 1, geo%nSpecies
        call asArray(lrN(iSp1),ctrl%dftbUInp%UJ(1:len(lrN(iSp1)),iSp1))
        allocate(iTmpN(len(liN(iSp1))))
        call asArray(liN(iSp1),iTmpN)
        ctrl%dftbUInp%niUJ(1:len(liN(iSp1)),iSp1) = iTmpN(:)
        deallocate(iTmpN)
        call destruct(lrN(iSp1))
        call destruct(liN(iSp1))
      end do
      allocate(ctrl%dftbUInp%iUJ(maxval(ctrl%dftbUInp%niUJ),&
          & maxval(ctrl%dftbUInp%nUJ),geo%nSpecies))
      ctrl%dftbUInp%iUJ(:,:,:) = 0
      do iSp1 = 1, geo%nSpecies
        do ii = 1, ctrl%dftbUInp%nUJ(iSp1)
          allocate(iTmpN(ctrl%dftbUInp%niUJ(ii,iSp1)))
          call get(li1N(iSp1),iTmpN,ii)
          ctrl%dftbUInp%iUJ(1:ctrl%dftbUInp%niUJ(ii,iSp1),ii,iSp1) = iTmpN(:)
          deallocate(iTmpN)
        end do
        call destruct(li1N(iSp1))
      end do

      deallocate(li1N)
      deallocate(lrN)
      deallocate(liN)

      ! check input values
      allocate(iTmpN(slako%orb%mShell))
      do iSp1 = 1, geo%nSpecies
        iTmpN = 0
        ! loop over number of blocks for that species
        do ii = 1, ctrl%dftbUInp%nUJ(iSp1)
          iTmpN(ctrl%dftbUInp%iUJ(1:ctrl%dftbUInp%niUJ(ii,iSp1),ii,iSp1)) = &
              & iTmpN(ctrl%dftbUInp%iUJ(1:ctrl%dftbUInp%niUJ(ii,iSp1),ii,iSp1)) + 1
        end do
        if (any(iTmpN(:)>1)) then
          write(stdout, *)'Multiple copies of shells present in OrbitalPotential!'
          write(stdout, "(A,A3,A,I2)") &
              & 'The count for the occurrence of shells of species ', &
              & trim(geo%speciesNames(iSp1)),' are:'
          write(stdout, *)iTmpN(1:slako%orb%nShell(iSp1))
          call abortProgram()
        end if
      end do
      deallocate(iTmpN)

    end if

    ! On-site
    call getChildValue(node, "OnSiteCorrection", value1, errStatus, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(value1)) then
      allocate(ctrl%onSiteElements(slako%orb%mShell, slako%orb%mShell, 2, geo%nSpecies))
      do iSp1 = 1, geo%nSpecies
        call getChildValue(child, trim(geo%speciesNames(iSp1))//"uu",&
            & ctrl%onSiteElements(:slako%orb%nShell(iSp1), :slako%orb%nShell(iSp1), 1, iSp1),&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child, trim(geo%speciesNames(iSp1))//"ud",&
            & ctrl%onSiteElements(:slako%orb%nShell(iSp1), :slako%orb%nShell(iSp1), 2, iSp1),&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end do
    end if

    ! Dispersion
    call getChildValue(node, "Dispersion", value1, errStatus, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(value1)) then
      allocate(ctrl%dispInp)
      call readDispersion(child, geo, ctrl%dispInp, ctrl%nrChrg, ctrl%tSCC, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! Solvation
    call getChildValue(node, "Solvation", value1, errStatus, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(value1)) then
      allocate(ctrl%solvInp)
      call readSolvation(child, geo, ctrl%solvInp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "RescaleSolvatedFields", ctrl%isSolvatedFieldRescaled, errStatus,&
          & .true.)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! Electronic constraints
    call getChildValue(node, "ElectronicConstraints", value1, errStatus, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true., list=.true.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(value1)) then
      allocate(ctrl%elecConstraintInp)
      call readElecConstraintInput(child, geo, ctrl%elecConstraintInp, ctrl%tSpin, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    if (ctrl%tLatOpt .and. .not. geo%tPeriodic) then
      @:RAISE_ERROR(errStatus, -1, "Lattice optimisation only applies for periodic structures.")
    end if

  #:if WITH_TRANSPORT
    call readElectrostatics(node, ctrl, geo, tp, poisson, errStatus)
    @:PROPAGATE_ERROR(errStatus)
  #:else
    call readElectrostatics(node, ctrl, geo, poisson, errStatus)
    @:PROPAGATE_ERROR(errStatus)
  #:endif

    ! Third order stuff
    ctrl%t3rd = .false.
    ctrl%t3rdFull = .false.
    if (ctrl%tSCC) then
      call getChildValue(node, "ThirdOrder", ctrl%t3rd, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(node, "ThirdOrderFull", ctrl%t3rdFull, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%t3rd .and. ctrl%t3rdFull) then
        call detailedError(node, "You must choose either ThirdOrder or ThirdOrderFull", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      if (ctrl%t3rd .and. ctrl%tShellResolved) then
        @:RAISE_ERROR(errStatus, -1, "Only full third-order DFTB is compatible with orbital&
            & resolved SCC")
      end if
      if (ctrl%t3rd .or. ctrl%t3rdFull) then
        call getChild(node, 'HubbardDerivs', child, errStatus, requested=.true.)
        @:PROPAGATE_ERROR(errStatus)
        allocate(ctrl%HubDerivs(slako%orb%mShell, geo%nSpecies))
        ctrl%hubDerivs(:,:) = 0.0_dp
        do iSp1 = 1, geo%nSpecies
          nShell = slako%orb%nShell(iSp1)
          if (ctrl%tShellResolved) then
            call getChildValue(child, geo%speciesNames(iSp1), ctrl%hubDerivs(1:nShell, iSp1),&
                & errStatus)
            @:PROPAGATE_ERROR(errStatus)
          else
            call getChildValue(child, geo%speciesNames(iSp1), ctrl%hubDerivs(1, iSp1), errStatus)
            @:PROPAGATE_ERROR(errStatus)
            ctrl%hubDerivs(2:nShell, iSp1) = ctrl%hubDerivs(1, iSp1)
          end if
        end do
        if (ctrl%t3rd) then
          allocate(ctrl%thirdOrderOn(geo%nAtom, 2))
          ctrl%thirdOrderOn(:,1) = 0.0_dp
          ctrl%thirdOrderOn(:,2) = ctrl%hubDerivs(1, geo%species)
        end if

        ! Halogen correction to the DFTB3 model
        block
          logical :: tHalogenInteraction
          integer :: iSp1, iSp2

          if (.not. geo%tPeriodic) then
            tHalogenInteraction = .false.
            iSp1Loop: do iSp1 = 1, geo%nSpecies
              if (any(geo%speciesNames(iSp1) == halogenXSpecies1)) then
                do iSp2 = 1, geo%nSpecies
                  if (any(geo%speciesNames(iSp2) == halogenXSpecies2)) then
                    tHalogenInteraction = .true.
                    exit iSp1Loop
                  end if
                end do
              end if
            end do iSp1Loop
            if (tHalogenInteraction) then
              call getChildValue(node, "HalogenXCorr", ctrl%tHalogenX, errStatus, .false.)
              @:PROPAGATE_ERROR(errStatus)
            end if
          end if
        end block

      end if
    end if

    call readDifferentiation(node, ctrl, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if (ctrl%tSCC) then
      ! Force type
      call readForceOptions(node, ctrl, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else
      ctrl%forceType = forceTypes%orig
    end if

    call readCustomisedHubbards(node, geo, slako%orb, ctrl%tShellResolved, ctrl%hubbU, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readDFTBHam


  !> Reads xTB-Hamiltonian
#:if WITH_TRANSPORT
  subroutine readXTBHam(node, ctrl, geo, tp, greendens, poisson, errStatus)
#:else
  subroutine readXTBHam(node, ctrl, geo, poisson, errStatus)
#:endif

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child, child2
    type(string) :: buffer, modifier
    type(string), allocatable :: searchPath(:)
    integer :: method, iSp1
    character(len=:), allocatable :: paramFile, paramTmp
    type(TOrbitals) :: orb

    ctrl%hamiltonian = hamiltonianTypes%xtb

    allocate(ctrl%tbliteInp)
    call ctrl%tbliteInp%setupGeometry(geo%nAtom, geo%species, geo%coords, geo%speciesNames,&
        & geo%latVecs)

    call getChild(node, "Method", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      call getChildValue(child, "", buffer, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      select case(unquote(char(buffer)))
      case default
        call detailedError(child, "Unknown method " // char(buffer) // " for xTB Hamiltonian",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      case("GFN1-xTB")
        method = tbliteMethod%gfn1xtb
      case("GFN2-xTB")
        method = tbliteMethod%gfn2xtb
      case("IPEA1-xTB")
        method = tbliteMethod%ipea1xtb
      end select
      call ctrl%tbliteInp%setupCalculator(method)
      ctrl%tbliteInp%info%name = trim(unquote(char(buffer)))
    else
      call getChildValue(node, "ParameterFile", value1, errStatus, "", child=child,&
          & allowEmptyValue=.true., dummyValue=.true.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(value1)) then
        call getChildValue(child, "", buffer, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        paramFile = trim(unquote(char(buffer)))
        call getParamSearchPath(searchPath)
        call findFile(searchPath, paramFile, paramTmp)
        if (allocated(paramTmp)) call move_alloc(paramTmp, paramFile)
        write(stdOut, '(a)') "Using parameter file '"//paramFile//"' for xTB Hamiltonian"
        call ctrl%tbliteInp%setupCalculator(paramFile)
      else
        call detailedError(node, "Either a Method or ParameterFile must be specified for xTB",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end if

    call getChildValue(node, "ShellResolvedSCC", ctrl%tShellResolved, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)

    ! SCC parameters
    call getChildValue(node, "SCC", ctrl%tSCC, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)
    ifSCC: if (ctrl%tSCC) then

      ! get charge mixing options etc.
      call readSccOptions(node, ctrl, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      !> TI-DFTB varibles for Delta DFTB
      call getChild(node, "NonAufbau", child, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child)) then
        ctrl%isNonAufbau = .true.
        call getChildValue(child, "SpinPurify", ctrl%isSpinPurify, errStatus, .true.)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child, "GroundGuess", ctrl%isGroundGuess, errStatus, .false.)
        @:PROPAGATE_ERROR(errStatus)
        ctrl%tSpin = .true.
        ctrl%t2Component = .false.
        ctrl%nrSpinPol = 0.0_dp
        ctrl%tSpinSharedEf = .false.
      else
        ctrl%isNonAufbau = .false.
      end if

    end if ifSCC

    ! Spin calculation
    if (ctrl%reksInp%reksAlg == reksTypes%noReks .and. .not.ctrl%isNonAufbau .and. ctrl%tSCC) then
    #:if WITH_TRANSPORT
      call readSpinPolarisation(node, ctrl, geo, tp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:else
      call readSpinPolarisation(node, ctrl, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:endif
    end if

    ! temporararily removed until debugged
    !if (.not. ctrl%tscc) then
    !  !! In a non-SCC calculation it is possible to upload charge shifts
    !  !! This is useful if the calculation can jump directly to the Analysis block
    !  call getChildValue(node, "ReadShifts", ctrl%tReadShifts, errStatus, .false.)
    !  @:PROPAGATE_ERROR(errStatus)
    !end if
    ctrl%tReadShifts = .false.

    ! External fields and potentials
    call readExternal(node, ctrl, geo, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Non-self-consistent spin-orbit coupling
    call ctrl%tbliteInp%setupOrbitals(geo%species, orb)
    call readSpinOrbit(node, ctrl, geo, orb, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Electronic solver
  #:if WITH_TRANSPORT
    call readSolver(node, ctrl, geo, tp, greendens, poisson, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if (tp%taskUpload) then
      ! Initialise, but unused
      ctrl%nrChrg =  0.0_dp
    else
      ! Charge
      call getChildValue(node, "Charge", ctrl%nrChrg, errStatus, 0.0_dp)
      @:PROPAGATE_ERROR(errStatus)
    end if
  #:else
    call readSolver(node, ctrl, geo, poisson, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Charge
    call getChildValue(node, "Charge", ctrl%nrChrg, errStatus, 0.0_dp)
    @:PROPAGATE_ERROR(errStatus)
  #:endif

    ! K-Points
    call readKPoints(node, ctrl, geo, ctrl%poorKSampling, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Dispersion
    call getChildValue(node, "Dispersion", value1, errStatus, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(value1)) then
      allocate(ctrl%dispInp)
      call readDispersion(child, geo, ctrl%dispInp, ctrl%nrChrg, ctrl%tSCC, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! Solvation
    call getChildValue(node, "Solvation", value1, errStatus, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(value1)) then
      allocate(ctrl%solvInp)
      call readSolvation(child, geo, ctrl%solvInp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "RescaleSolvatedFields", ctrl%isSolvatedFieldRescaled, errStatus,&
          & .true.)
      @:PROPAGATE_ERROR(errStatus)
    end if

    if (ctrl%tLatOpt .and. .not. geo%tPeriodic) then
      @:RAISE_ERROR(errStatus, -1, "Lattice optimisation only applies for periodic structures.")
    end if

  #:if WITH_TRANSPORT
    call readElectrostatics(node, ctrl, geo, tp, poisson, errStatus)
    @:PROPAGATE_ERROR(errStatus)
  #:else
    call readElectrostatics(node, ctrl, geo, poisson, errStatus)
    @:PROPAGATE_ERROR(errStatus)
  #:endif

    ! Third order stuff
    ctrl%t3rd = .true.
    ctrl%t3rdFull = .false.

    call readDifferentiation(node, ctrl, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    if (ctrl%tSCC) then
      ! Force type
      call readForceOptions(node, ctrl, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else
      ctrl%forceType = forceTypes%orig
    end if

  end subroutine readXTBHam


  !> Reads in settings for spin orbit enabled calculations
  subroutine readSpinOrbit(node, ctrl, geo, orb, errStatus)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Information about the orbitals of the species/atoms in the system
    class(TOrbitals), intent(in) :: orb

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child, child2
    type(string) :: modifier
    integer :: iSp

    call getChild(node, "SpinOrbit", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (.not. associated(child)) then
      ctrl%tSpinOrbit = .false.
      allocate(ctrl%xi(0,0))
    else
      if (ctrl%tSpin .and. .not. ctrl%t2Component) then
        @:RAISE_ERROR(errStatus, -1, "Spin-orbit coupling incompatible with collinear spin.")
      end if

      ctrl%tSpinOrbit = .true.
      ctrl%t2Component = .true.

      call getChildValue(child, "Dual", ctrl%tDualSpinOrbit, errStatus, .true.)
      @:PROPAGATE_ERROR(errStatus)

      allocate(ctrl%xi(orb%mShell,geo%nSpecies), source = 0.0_dp)
      do iSp = 1, geo%nSpecies
        call getChildValue(child, geo%speciesNames(iSp),&
            & ctrl%xi(:orb%nShell(iSp),iSp), errStatus, modifier=modifier, child=child2)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(modifier), energyUnits, child2, ctrl%xi(:orb%nShell(iSp),iSp),&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end do
    end if

  end subroutine readSpinOrbit


  !> Read in maximal angular momenta or selected shells
  subroutine readMaxAngularMomentum(node, geo, angShells, errStatus)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> List containing the angular momenta of the shells
    type(TListIntR1), allocatable, intent(out) :: angShells(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child, child2
    type(string) :: buffer
    integer :: iSp1, ii, jj, kk
    character(lc) :: strTmp
    integer :: nShell
    character(1) :: tmpCh
    type(TListString) :: lStr
    logical :: tShellIncl(4), tFound
    integer :: angShell(maxL+1), angShellOrdered(maxL+1)

    do ii = 1, maxL+1
      angShellOrdered(ii) = ii - 1
    end do
    call getChild(node, "MaxAngularMomentum", child, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    allocate(angShells(geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      call init(angShells(iSp1))
      call getChildValue(child, geo%speciesNames(iSp1), value1, errStatus, child=child2)
      @:PROPAGATE_ERROR(errStatus)
      call getNodeName(value1, buffer)
      select case(char(buffer))
      case("selectedshells")
        call init(lStr)
        call getChildValue(value1, "", lStr, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        do ii = 1, len(lStr)
          call get(lStr, strTmp, ii)
          strTmp = tolower(unquote(trim(strTmp)))
          if (len_trim(strTmp) > 4 .or. len_trim(strTmp) < 1) then
            call detailedError(value1, "Invalid shell selection '" // trim(strTmp)&
                & // "'. Nr. of selected shells must be between 1 and 4.", errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end if
          tShellIncl(:) = .false.
          nShell = len_trim(strTmp)
          do jj = 1, nShell
            tmpCh = strTmp(jj:jj)
            tFound = .false.
            do kk = 1, size(shellNames)
              if (tmpCh == trim(shellNames(kk))) then
                if (tShellIncl(kk)) then
                  call detailedError(value1, "Double selection of the same shell '" // tmpCh&
                      & // "' in shell selection block '" // trim(strTmp) // "'", errStatus)
                  @:PROPAGATE_ERROR(errStatus)
                end if
                tShellIncl(kk) = .true.
                angShell(jj) = kk - 1
                tFound = .true.
                exit
              end if
            end do
            if (.not. tFound) then
              call detailedError(value1, "Invalid shell name '" // tmpCh // "'", errStatus)
              @:PROPAGATE_ERROR(errStatus)
            end if
          end do
          call append(angShells(iSp1), angShell(1:nShell))
        end do
        call destruct(lStr)

      case(textNodeName)
        call getChildValue(child2, "", buffer, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        strTmp = unquote(char(buffer))
        do jj = 1, size(shellNames)
          if (trim(strTmp) == trim(shellNames(jj))) then
            call append(angShells(iSp1), angShellOrdered(:jj))
          end if
        end do
        if (len(angShells(iSp1)) < 1) then
          call detailedError(child2, "Invalid orbital name '" // trim(strTmp) // "'", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if

      case default
        call getNodeHSDName(value1, buffer)
        call detailedError(child2, "Invalid shell specification method '" // char(buffer) // "'",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end select
    end do

  end subroutine readMaxAngularMomentum


  !> Setup information about the orbitals of the species/atoms from angShell lists
  subroutine setupOrbitals(orb, geo, angShells, errStatus)

    !> Information about the orbitals of the species/atoms in the system
    class(TOrbitals), intent(out) :: orb

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> List containing the angular momenta of the shells,
    !> must be inout, since intoArray requires inout arguments
    type(TListIntR1), intent(inout) :: angShells(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: nShell, iSp1, iSh1, ii, jj, ind
    integer :: angShell(maxL+1)

    allocate(orb%nShell(geo%nSpecies))
    allocate(orb%nOrbSpecies(geo%nSpecies))
    allocate(orb%nOrbAtom(geo%nAtom))
    orb%mOrb = 0
    orb%mShell = 0
    do iSp1 = 1, geo%nSpecies
      orb%nShell(iSp1) = 0
      orb%nOrbSpecies(iSp1) = 0
      do ii = 1, len(angShells(iSp1))
        call intoArray(angShells(iSp1), angShell, nShell, ii)
        orb%nShell(iSp1) = orb%nShell(iSp1) + nShell
        do jj = 1, nShell
          orb%nOrbSpecies(iSp1) = orb%nOrbSpecies(iSp1) &
              &+ 2 * angShell(jj) + 1
        end do
      end do
    end do
    orb%mShell = maxval(orb%nShell)
    orb%mOrb = maxval(orb%nOrbSpecies)
    orb%nOrbAtom(:) = orb%nOrbSpecies(geo%species(:))
    orb%nOrb = sum(orb%nOrbAtom)

    allocate(orb%angShell(orb%mShell, geo%nSpecies))
    allocate(orb%iShellOrb(orb%mOrb, geo%nSpecies))
    allocate(orb%posShell(orb%mShell+1, geo%nSpecies))
    orb%angShell(:,:) = 0
    do iSp1 = 1, geo%nSpecies
      ind = 1
      iSh1 = 1
      do ii = 1, len(angShells(iSp1))
        call intoArray(angShells(iSp1), angShell, nShell, ii)
        do jj = 1, nShell
          orb%posShell(iSh1, iSp1) = ind
          orb%angShell(iSh1, iSp1) = angShell(jj)
          orb%iShellOrb(ind:ind+2*angShell(jj), iSp1) = iSh1
          ind = ind + 2 * angShell(jj) + 1
          iSh1 = iSh1 + 1
        end do
        orb%posShell(iSh1, iSp1) = ind
      end do
    end do

  end subroutine setupOrbitals


#:if WITH_TRANSPORT
  subroutine readElectrostatics(node, ctrl, geo, tp, poisson, errStatus)
#:else
  subroutine readElectrostatics(node, ctrl, geo, poisson, errStatus)
#:endif

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp
  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child
    type(string) :: buffer

    ctrl%tPoisson = .false.

    ! Read in which kind of electrostatics method to use.
    call getChildValue(node, "Electrostatics", value1, errStatus, "GammaFunctional", child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)

    select case (char(buffer))

    case ("gammafunctional")
    #:if WITH_TRANSPORT
      if (tp%taskUpload .and. ctrl%tSCC) then
        call detailedError(value1, "GammaFunctional not available, if you upload contacts in an SCC&
            & calculation.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    #:endif

    case ("poisson")
      if (.not. withPoisson) then
        call detailedError(value1, "Poisson not available as binary was built without the Poisson&
            &-solver", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        ctrl%tPoisson = .true.
        #:if WITH_TRANSPORT
        call readPoisson(value1, poisson, geo%tPeriodic, tp, geo%latVecs, ctrl%updateSccAfterDiag,&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
        #:else
        call readPoisson(value1, poisson, geo%tPeriodic, geo%latVecs, ctrl%updateSccAfterDiag,&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
        #:endif
      #:endblock

    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Unknown electrostatics '" // char(buffer) // "'", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readElectrostatics


  !> Spin calculation
#:if WITH_TRANSPORT
  subroutine readSpinPolarisation(node, ctrl, geo, tp, errStatus)
#:else
  subroutine readSpinPolarisation(node, ctrl, geo, errStatus)
#:endif

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp
  #:endif

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child
    type(string) :: buffer

    call renameChildren(node, "SpinPolarization", "SpinPolarisation")
    call getChildValue(node, "SpinPolarisation", value1, errStatus, "", child=child,&
        & allowEmptyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName2(value1, buffer)
    select case(char(buffer))
    case ("")
      ctrl%tSpin = .false.
      ctrl%t2Component = .false.
      ctrl%nrSpinPol = 0.0_dp

    case ("colinear", "collinear")
      ctrl%tSpin = .true.
      ctrl%t2Component = .false.
      call getChildValue(value1, 'UnpairedElectrons', ctrl%nrSpinPol, errStatus, 0.0_dp)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, 'RelaxTotalSpin', ctrl%tSpinSharedEf, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
      if (.not. ctrl%tReadChrg) then
        call getInitialSpins(value1, geo, 1, ctrl%initialSpins, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

    case ("noncolinear", "noncollinear")
      ctrl%tSpin = .true.
      ctrl%t2Component = .true.
      if (.not. ctrl%tReadChrg) then
        call getInitialSpins(value1, geo, 3, ctrl%initialSpins, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Invalid spin polarisation type '" // char(buffer) // "'",&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readSpinPolarisation


  ! External field(s) and potential(s)
  subroutine readExternal(node, ctrl, geo, errStatus)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child, child2, child3
    type(fnodeList), pointer :: children
    type(string) :: modifier, buffer, buffer2
    real(dp) :: rTmp
    type(TFileDescr) :: file
    integer :: ind, ii, iErr, nElem
    real(dp), allocatable :: tmpR1(:), tmpR2(:,:)
    type(TListRealR2) :: lCharges
    type(TListRealR1) :: lBlurs, lr1
    type(TListReal) :: lr
    type(TListInt) :: li

    call getChildValue(node, "ElectricField", value1, errStatus, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true., list=.true.)
    @:PROPAGATE_ERROR(errStatus)

    ! external applied field
    call getChild(child, "External", child2, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child2)) then
      allocate(ctrl%electricField)
      ctrl%tMulliken = .true.
      call getChildValue(child2, "Strength", ctrl%electricField%EFieldStrength, errStatus,&
          & modifier=modifier, child=child3)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), EFieldUnits, child3, ctrl%electricField%EFieldStrength,&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child2, "Direction", ctrl%electricField%EfieldVector, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (sum(ctrl%electricField%EfieldVector**2) < 1e-8_dp) then
        call detailedError(child2, "Vector too small", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      else
        ctrl%electricField%EfieldVector = ctrl%electricField%EfieldVector&
            & / sqrt(sum(ctrl%electricField%EfieldVector**2))
      end if
      call getChildValue(child2, "Frequency", ctrl%electricField%EFieldOmega, errStatus, 0.0_dp,&
          & modifier=modifier, child=child3)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), freqUnits, child3, ctrl%electricField%EFieldOmega,&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%electricField%EFieldOmega > 0.0) then
        ! angular frequency
        ctrl%electricField%EFieldOmega = 2.0_dp * pi * ctrl%electricField%EFieldOmega
        ctrl%electricField%isTDEfield = .true.
      else
        ctrl%electricField%isTDEfield = .false.
        ctrl%electricField%EFieldOmega = 0.0_dp
      end if
      ctrl%electricField%EfieldPhase = 0
      if (ctrl%electricField%isTDEfield) then
        call getChildValue(child2, "Phase", ctrl%electricField%EfieldPhase, errStatus, 0)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end if

    ctrl%nExtChrg = 0
    if (ctrl%hamiltonian == hamiltonianTypes%dftb) then

      call getChildren(child, "PointCharges", children)
      if (getLength(children) > 0) then
        ! Point charges present
        if (.not.ctrl%tSCC) then
          @:RAISE_ERROR(errStatus, -1, "External charges can only be used in an SCC calculation")
        end if
        call init(lCharges)
        call init(lBlurs)
        ctrl%nExtChrg = 0
        do ii = 1, getLength(children)
          call getItem1(children, ii, child2)
          call getChildValue(child2, "CoordsAndCharges", value1, errStatus, modifier=modifier,&
              & child=child3)
          @:PROPAGATE_ERROR(errStatus)
          call getNodeName(value1, buffer)
          select case(char(buffer))
          case (textNodeName)
            call init(lr1)
            call getChildValue(child3, "", 4, lr1, errStatus, modifier=modifier)
            @:PROPAGATE_ERROR(errStatus)
            allocate(tmpR2(4, len(lr1)))
            call asArray(lr1, tmpR2)
            ctrl%nExtChrg = ctrl%nExtChrg + len(lr1)
            call destruct(lr1)
          case ("directread")
            call getChildValue(value1, "Records", ind, errStatus)
            @:PROPAGATE_ERROR(errStatus)
            call getChildValue(value1, "File", buffer2, errStatus)
            @:PROPAGATE_ERROR(errStatus)
            allocate(tmpR2(4, ind))
            call openFile(file, unquote(char(buffer2)), mode="r", iostat=iErr)
            if (iErr /= 0) then
              call detailedError(value1, "Could not open file '" // trim(unquote(char(buffer2)))&
                  & // "' for direct reading", errStatus)
              @:PROPAGATE_ERROR(errStatus)
            end if
            read(file%unit, *, iostat=iErr) tmpR2
            if (iErr /= 0) then
              call detailedError(value1, "Error during direct reading '"&
                  & // trim(unquote(char(buffer2))) // "'", errStatus)
              @:PROPAGATE_ERROR(errStatus)
            end if
            call closeFile(file)
            ctrl%nExtChrg = ctrl%nExtChrg + ind
          case default
            call detailedError(value1, "Invalid block name", errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end select
          call convertUnitHsd(char(modifier), lengthUnits, child3, tmpR2(1:3,:), errStatus)
          @:PROPAGATE_ERROR(errStatus)
          call append(lCharges, tmpR2)
          call getChildValue(child2, "GaussianBlurWidth", rTmp, errStatus, 0.0_dp,&
              & modifier=modifier, child=child3)
          @:PROPAGATE_ERROR(errStatus)
          if (rTmp < 0.0_dp) then
            call detailedError(child3, "Gaussian blur width may not be negative", errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end if
          call convertUnitHsd(char(modifier), lengthUnits, child3, rTmp, errStatus)
          @:PROPAGATE_ERROR(errStatus)
          allocate(tmpR1(size(tmpR2, dim=2)))
          tmpR1(:) = rTmp
          call append(lBlurs, tmpR1)
          deallocate(tmpR1)
          deallocate(tmpR2)
        end do

        allocate(ctrl%extChrg(4, ctrl%nExtChrg))
        ind = 1
        do ii = 1, len(lCharges)
          call intoArray(lCharges, ctrl%extChrg(:, ind:), nElem, ii)
          ind = ind + nElem
        end do
        call destruct(lCharges)

        allocate(ctrl%extChrgBlurWidth(ctrl%nExtChrg))
        ind = 1
        do ii = 1, len(lBlurs)
          call intoArray(lBlurs, ctrl%extChrgBlurWidth(ind:), nElem, ii)
          ind = ind + nElem
        end do
        call destruct(lBlurs)
        call destroyNodeList(children)
      end if

    else

      call getChildren(child, "PointCharges", children)
      if (getLength(children) > 0) then
        call detailedError(child, "External charges are not currently supported for this model",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

    end if

    call getChild(node, "AtomSitePotential", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      allocate(ctrl%atomicExtPotential)

      call getChild(child, "Net", child2, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child2)) then
        ! onsites
        ctrl%tNetAtomCharges = .true.
        call init(li)
        call init(lr)
        call getChildValue(child2, "Atoms", li, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child2, "Vext", lr, errStatus, modifier=modifier, child=child3)
        @:PROPAGATE_ERROR(errStatus)
        if (len(li) /= len(lr)) then
          call detailedError(child2, "Mismatch in number of sites and potentials", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        allocate(ctrl%atomicExtPotential%iAtOnSite(len(li)))
        call asArray(li, ctrl%atomicExtPotential%iAtOnSite)
        allocate(ctrl%atomicExtPotential%VextOnSite(len(lr)))
        call asArray(lr, ctrl%atomicExtPotential%VextOnSite)
        call convertUnitHsd(char(modifier), energyUnits, child3,&
            & ctrl%atomicExtPotential%VextOnSite, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call destruct(li)
        call destruct(lr)
      end if

      call getChild(child, "Gross", child2, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child2)) then
        ! atomic
        call init(li)
        call init(lr)
        call getChildValue(child2, "Atoms", li, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child2, "Vext", lr, errStatus, modifier=modifier, child=child3)
        @:PROPAGATE_ERROR(errStatus)
        if (len(li) /= len(lr)) then
          call detailedError(child2, "Mismatch in number of sites and potentials", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        allocate(ctrl%atomicExtPotential%iAt(len(li)))
        call asArray(li, ctrl%atomicExtPotential%iAt)
        allocate(ctrl%atomicExtPotential%Vext(len(lr)))
        call asArray(lr, ctrl%atomicExtPotential%Vext)
        call convertUnitHsd(char(modifier), energyUnits, child3, ctrl%atomicExtPotential%Vext,&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call destruct(li)
        call destruct(lr)
      end if

      if (.not.allocated(ctrl%atomicExtPotential%iAt)&
          & .and. .not.allocated(ctrl%atomicExtPotential%iAtOnSite)) then
        call detailedError(child, "No atomic potentials specified", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

    end if

  end subroutine readExternal


  !> Filling of electronic levels
  subroutine readFilling(node, ctrl, geo, temperatureDefault, errStatus)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to test for periodicity
    type(TGeometry), intent(in) :: geo

    !> Default temperature for filling
    real(dp), intent(in) :: temperatureDefault

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child, child2, child3, field
    type(string) :: buffer, modifier
    character(lc) :: errorStr

    call getChildValue(node, "Filling", value1, errStatus, "Fermi", child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)

    select case (char(buffer))
    case ("fermi")
      ctrl%iDistribFn = fillingTypes%Fermi ! Fermi function
    case ("methfesselpaxton")
      ! Set the order of the Methfessel-Paxton step function approximation, defaulting to 2nd order
      call getChildValue(value1, "Order", ctrl%iDistribFn, errStatus, 2)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%iDistribFn < 1) then
        call getNodeHSDName(value1, buffer)
        write(errorStr, "(A,A,A,I4)")"Unsuported filling mode for '", &
            & char(buffer),"' :",ctrl%iDistribFn
        call detailedError(child, errorStr, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      ctrl%iDistribFn = fillingTypes%Methfessel + ctrl%iDistribFn
    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Invalid filling method '" //char(buffer)// "'", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    if (.not. ctrl%tSetFillingTemp) then
      call getChildValue(value1, "Temperature", ctrl%tempElec, errStatus, temperatureDefault,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, field, ctrl%tempElec, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%tempElec < minTemp) then
        ctrl%tempElec = minTemp
      end if
    end if

    call getChild(value1, "FixedFermiLevel", child2, errStatus, modifier=modifier,&
        & requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    ctrl%tFixEf = associated(child2)
    if (ctrl%tFixEf) then
      if (ctrl%tSpin .and. .not.ctrl%t2Component) then
        allocate(ctrl%Ef(2))
      else
        allocate(ctrl%Ef(1))
      end if
      call getChildValue(child2, "", ctrl%Ef, errStatus, modifier=modifier, child=child3)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, child3, ctrl%Ef, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    if (geo%tPeriodic .and. .not.ctrl%tFixEf) then
      call getChildValue(value1, "IndependentKFilling", ctrl%tFillKSep, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine readFilling


  !> Electronic Solver
#:if WITH_TRANSPORT
  subroutine readSolver(node, ctrl, geo, tp, greendens, poisson, errStatus)
#:else
  subroutine readSolver(node, ctrl, geo, poisson, errStatus)
#:endif

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

  #:endif

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child
    type(string) :: buffer, modifier

    integer :: iTmp

    ! Electronic solver
    call getChildValue(node, "Solver", value1, errStatus, "RelativelyRobust")
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)

    select case(char(buffer))

    case ("qr")
      ctrl%solver%isolver = electronicSolverTypes%qr

    case ("divideandconquer")
      ctrl%solver%isolver = electronicSolverTypes%divideandconquer

    case ("relativelyrobust")
      ctrl%solver%isolver = electronicSolverTypes%relativelyrobust

  #:if WITH_MAGMA
    case ("magma")
      ctrl%solver%isolver = electronicSolverTypes%magma_gvd
  #:endif

    case ("elpa")
      allocate(ctrl%solver%elsi)
      call getChildValue(value1, "Sparse", ctrl%solver%elsi%elsiCsr, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%solver%elsi%elsiCsr) then
        ctrl%solver%isolver = electronicSolverTypes%elpadm
      else
        ctrl%solver%isolver = electronicSolverTypes%elpa
      end if
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      call getChildValue(value1, "Mode", ctrl%solver%elsi%elpaSolver, errStatus, 2)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "Autotune", ctrl%solver%elsi%elpaAutotune, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "Gpu", ctrl%solver%elsi%elpaGpu, errStatus, .false., child=child)
      @:PROPAGATE_ERROR(errStatus)
      #:if not WITH_GPU
        if (ctrl%solver%elsi%elpaGpu) then
          call detailedError(child, "DFTB+ must be compiled with GPU support in order to enable&
              & the GPU acceleration for the ELPA solver", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      #:endif

    case ("omm")
      ctrl%solver%isolver = electronicSolverTypes%omm
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      call getChildValue(value1, "nIterationsELPA", ctrl%solver%elsi%ommIterationsElpa, errStatus,&
          & 5)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "Tolerance", ctrl%solver%elsi%ommTolerance, errStatus, 1.0E-10_dp)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "Choleskii", ctrl%solver%elsi%ommCholesky, errStatus, .true.)
      @:PROPAGATE_ERROR(errStatus)

    case ("pexsi")
      ctrl%solver%isolver = electronicSolverTypes%pexsi
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
    #:if ELSI_VERSION > 2.5
      call getChildValue(value1, "Method", ctrl%solver%elsi%pexsiMethod, errStatus, 3)
      @:PROPAGATE_ERROR(errStatus)
    #:else
      call getChildValue(value1, "Method", ctrl%solver%elsi%pexsiMethod, errStatus, 2)
      @:PROPAGATE_ERROR(errStatus)
    #:endif
      select case(ctrl%solver%elsi%pexsiMethod)
      case(1)
        iTmp = 60
      case(2)
        iTmp = 20
      case(3)
        iTmp = 30
      end select
      call getChildValue(value1, "Poles", ctrl%solver%elsi%pexsiNPole, errStatus, iTmp)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%solver%elsi%pexsiNPole < 10) then
        call detailedError(value1, "Too few PEXSI poles", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      select case(ctrl%solver%elsi%pexsiMethod)
      case(1)
        if (mod(ctrl%solver%elsi%pexsiNPole,10) /= 0 .or. ctrl%solver%elsi%pexsiNPole > 120) then
          call detailedError(value1, "Unsupported number of PEXSI poles for method 1", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      case(2,3)
        if (mod(ctrl%solver%elsi%pexsiNPole,5) /= 0 .or. ctrl%solver%elsi%pexsiNPole > 40) then
          call detailedError(value1, "Unsupported number of PEXSI poles for this method", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      end select
      call getChildValue(value1, "ProcsPerPole", ctrl%solver%elsi%pexsiNpPerPole, errStatus, 1)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "muPoints", ctrl%solver%elsi%pexsiNMu, errStatus, 2)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "SymbolicFactorProcs", ctrl%solver%elsi%pexsiNpSymbo, errStatus, 1)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "SpectralRadius", ctrl%solver%elsi%pexsiDeltaE, errStatus,&
          & 10.0_dp, modifier=modifier, child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, child, ctrl%solver%elsi%pexsiDeltaE,&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case ("ntpoly")
      ctrl%solver%isolver = electronicSolverTypes%ntpoly
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      if (ctrl%tSpin) then
        call detailedError(value1, "Solver does not currently support spin polarisation", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getChildValue(value1, "PurificationMethod", ctrl%solver%elsi%ntpolyMethod, errStatus, 2)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "Tolerance", ctrl%solver%elsi%ntpolyTolerance, errStatus,&
          & 1.0E-5_dp)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "Truncation", ctrl%solver%elsi%ntpolyTruncation, errStatus,&
          & 1.0E-10_dp)
      @:PROPAGATE_ERROR(errStatus)

  #:if WITH_TRANSPORT
    case ("greensfunction")
      ctrl%solver%isolver = electronicSolverTypes%GF
      ! need electronic temperature to be read for this solver:
      call readElectronicFilling(node, ctrl, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (tp%defined .and. .not.tp%taskUpload) then
        call detailederror(node, "greensfunction solver cannot be used when task&
            & = contactHamiltonian", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call readGreensFunction(value1, greendens, tp, ctrl%tempElec, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      ! fixEf also avoids checks of total charge in initQFromFile
      ctrl%tFixEf = .true.
    case ("transportonly")
      if (tp%defined .and. .not.tp%taskUpload) then
        call detailederror(node, "transportonly cannot be used when task = contactHamiltonian",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call readGreensFunction(value1, greendens, tp, ctrl%tempElec, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      ctrl%solver%isolver = electronicSolverTypes%OnlyTransport
      ctrl%tFixEf = .true.
  #:endif

    case default
      call detailedError(value1, "Unknown electronic solver", errStatus)
      @:PROPAGATE_ERROR(errStatus)

    end select

    if ((ctrl%solver%isolver == electronicSolverTypes%omm .or.&
        & ctrl%solver%isolver == electronicSolverTypes%pexsi ) .and. .not.ctrl%tSpinSharedEf&
        & .and. ctrl%tSpin .and. .not. ctrl%t2Component) then
      call detailedError(value1, "This solver currently requires spin values to be relaxed",&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    if (ctrl%solver%isolver == electronicSolverTypes%pexsi .and. .not.withPEXSI) then
      @:RAISE_ERROR(errStatus, -1, "Not compiled with PEXSI support via ELSI")
    end if
    if (any(ctrl%solver%isolver == [electronicSolverTypes%elpa, electronicSolverTypes%omm,&
        & electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly])) then
      if (.not.withELSI) then
        @:RAISE_ERROR(errStatus, -1, "Not compiled with ELSI supported solvers")
      end if
    end if

    if (any(ctrl%solver%isolver == [electronicSolverTypes%omm, electronicSolverTypes%pexsi,&
        & electronicSolverTypes%ntpoly])) then
      call getChildValue(value1, "Sparse", ctrl%solver%elsi%elsiCsr, errStatus, .true.)
      @:PROPAGATE_ERROR(errStatus)
      if (.not.ctrl%solver%elsi%elsiCsr) then
        if (any(ctrl%solver%isolver == [electronicSolverTypes%pexsi,electronicSolverTypes%ntpoly]))&
            & then
          call getChildValue(value1, "Threshold", ctrl%solver%elsi%elsi_zero_def, errStatus,&
              & 1.0E-15_dp)
          @:PROPAGATE_ERROR(errStatus)
        end if
      end if
    end if

  #:if WITH_TRANSPORT
    if (all(ctrl%solver%isolver /= [electronicSolverTypes%GF,electronicSolverTypes%OnlyTransport])&
        & .and. tp%taskUpload) then
      call detailedError(value1, "Eigensolver incompatible with transport calculation&
          & (GreensFunction or TransportOnly required)", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
  #:endif
  end subroutine readSolver


  !> K-Points
  subroutine readKPoints(node, ctrl, geo, poorKSampling, errStatus)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    !> Is this k-point grid usable to integrate properties like the energy, charges, ...?
    logical, intent(out) :: poorKSampling

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: ii
    character(lc) :: errorStr

    ! Assume SCC can has usual default number of steps if needed
    poorKSampling = .false.

    ! K-Points
    if (geo%tPeriodic) then
      call getEuclideanKSampling(poorKSampling, ctrl, node, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else if (geo%tHelical) then
      call getHelicalKSampling(poorKSampling, ctrl, node, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call maxSelfConsIterations(node, ctrl, "MaxSCCIterations", ctrl%maxSccIter, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    ! Eventually, perturbation routines should also have restart reads:
    if (ctrl%poorKSampling .and. ctrl%tSCC .and. .not.ctrl%tReadChrg) then
      call warning("It is strongly suggested you use the ReadInitialCharges option.")
    end if

  end subroutine readKPoints


  !> Set the maximum number of SCC cycles, depending on k-point behaviour
  subroutine maxSelfConsIterations(node, ctrl, label, maxSccIter, errStatus)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Name of the tag
    character(*), intent(in) :: label

    !> Number of self-consistent iterations
    integer, intent(out) :: maxSccIter

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    ! string for error return
    character(lc) :: warningStr

    integer :: ii

    maxSccIter = 1
    if (ctrl%tSCC) then
      if (ctrl%poorKSampling) then
        ! prevent full SCC with these points
        ii = 1
      else
        ii = 100
      end if
      call getChildValue(node, trim(label), maxSccIter, errStatus, ii)
      @:PROPAGATE_ERROR(errStatus)
    end if

    if (ctrl%poorKSampling .and. maxSccIter /= 1) then
      write(warningStr, "(A,I3)") "A self-consistent cycle with these k-points probably will&
          & not correctly calculate many properties, maximum iterations set to:", maxSccIter
      call warning(warningStr)
    end if

  end subroutine maxSelfConsIterations


  !> K-points in Euclidean space
  subroutine getEuclideanKSampling(poorKSampling, ctrl, node, geo, errStatus)

    !> Is this k-point grid usable to integrate properties like the energy, charges, ...?
    logical, intent(out) :: poorKSampling

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child
    type(string) :: buffer, modifier
    integer :: ind, ii, jj, kk
    real(dp) :: coeffsAndShifts(3, 4)
    real(dp) :: rTmp3(3)
    type(TListIntR1) :: li1
    type(TListRealR1) :: lr1
    integer, allocatable :: tmpI1(:)
    real(dp), allocatable :: kpts(:,:)

    call getChildValue(node, "KPointsAndWeights", value1, errStatus, child=child,&
        & modifier=modifier)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    select case(char(buffer))

    case ("supercellfolding")
      poorKSampling = .false.
      if (len(modifier) > 0) then
        call detailedError(child, "No modifier is allowed, if the SupercellFolding scheme is&
            & used.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getChildValue(value1, "", coeffsAndShifts, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (abs(determinant33(coeffsAndShifts(:,1:3))) - 1.0_dp < -1e-6_dp) then
        call detailedError(value1, "Determinant of the supercell matrix must be greater than 1",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      if (any(abs(modulo(coeffsAndShifts(:,1:3) + 0.5_dp, 1.0_dp) - 0.5_dp) &
          &> 1e-6_dp)) then
        call detailedError(value1, "The components of the supercell matrix must be integers.",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      if (.not.ctrl%tSpinOrbit) then
        call getSuperSampling(coeffsAndShifts(:,1:3), modulo(coeffsAndShifts(:,4), 1.0_dp),&
            & ctrl%kPoint, ctrl%kWeight, reduceByInversion=.true.)
      else
        call getSuperSampling(coeffsAndShifts(:,1:3), modulo(coeffsAndShifts(:,4), 1.0_dp),&
            & ctrl%kPoint, ctrl%kWeight, reduceByInversion=.false.)
      end if
      ctrl%nKPoint = size(ctrl%kPoint, dim=2)

    case ("klines")
      ! probably unable to integrate charge for SCC
      poorKSampling = .true.
      call init(li1)
      call init(lr1)
      call getChildValue(value1, "", 1, li1, 3, lr1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (len(li1) < 1) then
        call detailedError(value1, "At least one line must be specified.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      allocate(tmpI1(len(li1)))
      allocate(kpts(3, 0:len(lr1)))
      call asVector(li1, tmpI1)
      call asArray(lr1, kpts(:,1:len(lr1)))
      kpts(:,0) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
      call destruct(li1)
      call destruct(lr1)
      if (any(tmpI1 < 0)) then
        call detailedError(value1, "Interval steps must be greater equal to zero.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      ctrl%nKPoint = sum(tmpI1)
      if (ctrl%nKPoint < 1) then
        call detailedError(value1, "Sum of the interval steps must be greater than zero.",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      ii = 1
      do while (tmpI1(ii) == 0)
        ii = ii + 1
      end do
      allocate(ctrl%kPoint(3, ctrl%nKPoint))
      allocate(ctrl%kWeight(ctrl%nKPoint))
      ind = 1
      do jj = ii, size(tmpI1)
        if (tmpI1(jj) == 0) then
          cycle
        end if
        rTmp3(:) = (kpts(:,jj) - kpts(:,jj-1)) / real(tmpI1(jj), dp)
        do kk = 1, tmpI1(jj)
          ctrl%kPoint(:,ind) = kpts(:,jj-1) + real(kk, dp) * rTmp3
          ind = ind + 1
        end do
      end do
      ctrl%kWeight(:) = 1.0_dp
      if (len(modifier) > 0) then
        select case (tolower(char(modifier)))
        case ("relative")
        case ("absolute")
          ctrl%kPoint(:,:) =  matmul(transpose(geo%latVecs), ctrl%kPoint)
          kpts(:,:) = matmul(transpose(geo%latVecs), kpts)
        case default
          call detailedError(child, "Invalid modifier: '" // char(modifier) // "'", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end select
      end if
      deallocate(tmpI1)
      deallocate(kpts)

    case (textNodeName)

      ! no idea, but assume user knows what they are doing
      poorKSampling = .false.

      call init(lr1)
      call getChildValue(child, "", 4, lr1, errStatus, modifier=modifier)
      @:PROPAGATE_ERROR(errStatus)
      if (len(lr1) < 1) then
        call detailedError(child, "At least one k-point must be defined.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      ctrl%nKPoint = len(lr1)
      allocate(kpts(4, ctrl%nKPoint))
      call asArray(lr1, kpts)
      call destruct(lr1)
      if (len(modifier) > 0) then
        select case (tolower(char(modifier)))
        case ("relative")
          continue
        case ("absolute")
          kpts(1:3,:) =  matmul(transpose(geo%latVecs), kpts(1:3,:))
        case default
          call detailedError(child, "Invalid modifier: '" // char(modifier) // "'", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end select
      end if
      allocate(ctrl%kPoint(3, ctrl%nKPoint))
      allocate(ctrl%kWeight(ctrl%nKPoint))
      ctrl%kPoint(:,:) = kpts(1:3, :)
      ctrl%kWeight(:) = kpts(4, :)
      deallocate(kpts)
    case default
      call detailedError(value1, "Invalid K-point scheme", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine getEuclideanKSampling


  !> K-points for helical boundaries
  subroutine getHelicalKSampling(poorKSampling, ctrl, node, geo, errStatus)

    !> Is this k-point grid usable to integrate properties like the energy, charges, ...?
    logical, intent(out) :: poorKSampling

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(string) :: buffer
    type(fnode), pointer :: value1, child
    type(TListRealR1) :: lr1
    real(dp):: rTmp3(3), rTmp22(2,2)
    integer :: iTmp, iTmp2(2), kk, ii, jj
    real(dp), allocatable :: kPts(:,:)
    character(lc) :: errorStr

    ! assume the user knows what they are doing
    poorKSampling = .false.

    call getChildValue(node, "KPointsAndWeights", value1, errStatus, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case ("helicaluniform")
      call getChildValue(value1, "", rTmp3(:2), errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (abs(modulo(rTmp3(1) + 0.5_dp, 1.0_dp) - 0.5_dp) > 1e-6_dp) then
        call detailedError(value1, "The k-point grid must be integer values.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      iTmp = nint(rTmp3(1))
      if (iTmp < 1) then
        call detailedError(node, "Number of grid points must be above 0", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      if (.not.ctrl%tSpinOrbit) then
        ctrl%nKPoint = iTmp * nint(geo%latvecs(3,1))
        allocate(ctrl%kPoint(2, ctrl%nKPoint))
        ctrl%kPoint(:,:) = 0.0_dp
        allocate(ctrl%kWeight(ctrl%nKPoint))
        ctrl%kWeight(:) = 1.0_dp / real(iTmp,dp)
        do ii = 0, iTmp-1
          ctrl%kPoint(1,ii+1) = ii * 0.5_dp*ctrl%kWeight(ii+1) + 0.5_dp*rTmp3(2)/rTmp3(1)
        end do
        ctrl%kWeight(:) = 1.0_dp / real(ctrl%nKPoint,dp)
        do ii = 2, nint(geo%latvecs(3,1))
          ctrl%kPoint(1,(ii-1)*iTmp+1:ii*iTmp) = ctrl%kPoint(1,1:iTmp)
          ctrl%kPoint(2,(ii-1)*iTmp+1:ii*iTmp) = real(ii-1,dp)/nint(geo%latvecs(3,1))
        end do
      else
        @:RAISE_ERROR(errStatus, -1, "Helical boundaries not yet added for spin-orbit")
      end if
    case ("helicalsampled")
      call getChildValue(value1, "", rTmp22, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      iTmp2 = nint(rTmp22(:,1))
      if (any(abs(iTmp2-rTmp22(:,1)) > 1e-6_dp)) then
        call detailedError(value1, "The k-point grid must be integers.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      if (any(iTmp2 < 1)) then
        call detailedError(node, "Number of grid points must be above 0", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      if (iTmp2(2) > nint(geo%latvecs(3,1))) then
        write(errorStr, '("The k-point grid for the helix rotational operation (",I0,&
            & ") is larger than the rotation order (C_",I0,").")') iTmp2(2), nint(geo%latvecs(3,1))
        call detailedError(node, errorStr, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      if (mod(nint(geo%latvecs(3,1)),iTmp2(2)) /= 0) then
        write(errorStr, '("The k-point grid for the helix rotational operation (n_k=",I0,&
            & ") is not a divisor of the rotation order (C_",I0,").")') iTmp2(2),&
            & nint(geo%latvecs(3,1))
        call detailedError(node, errorStr, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      if (abs(rTmp22(2,2) * nint(geo%latvecs(3,1)) - nint(rTmp22(2,2) * nint(geo%latvecs(3,1))))&
          & > epsilon(1.0_dp)) then
        write(errorStr, '("The shift of the k-points along the rotation is incommensurate, it must&
            & be an integer multiple of 1/",I0)') nint(geo%latvecs(3,1))
        call detailedError(node, errorStr, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      if (.not.ctrl%tSpinOrbit) then
        ctrl%nKPoint = product(iTmp2)
        allocate(ctrl%kPoint(2, ctrl%nKPoint))
        ctrl%kPoint(:,:) = 0.0_dp
        allocate(ctrl%kWeight(ctrl%nKPoint))

        kk = 1
        do ii = 0, iTmp2(1)-1
          do jj = 0, iTmp2(2)-1
            ctrl%kPoint(1,kk) = ii * 0.5_dp / rTmp22(1,1) + 0.5_dp*rTmp22(1,2)/rTmp22(1,1)
            ctrl%kPoint(2,kk) = mod(jj * 1.0_dp / rTmp22(2,1) + rTmp22(2,2), 1.0_dp)
            kk = kk + 1
          end do
        end do

        ctrl%kWeight(:) = 1.0_dp / real(ctrl%nKPoint,dp)

      else
        @:RAISE_ERROR(errStatus, -1, "Helical boundaries not yet added for spin-orbit")
      end if

    case (textNodeName)

      call init(lr1)
      call getChildValue(child, "", 3, lr1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (len(lr1) < 1) then
        call detailedError(child, "At least one k-point must be defined.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      ctrl%nKPoint = len(lr1)
      allocate(kpts(3, ctrl%nKPoint))
      call asArray(lr1, kpts)
      call destruct(lr1)
      allocate(ctrl%kPoint(2, ctrl%nKPoint))
      allocate(ctrl%kWeight(ctrl%nKPoint))
      ! first two are point values
      ctrl%kPoint(:2,:) = kpts(:2, :)
      ! test if the second k-point is commensurate with the C_n operation
      if (any(abs(kpts(2,:)*nint(geo%latvecs(3,1)) - nint(kpts(2,:) * nint(geo%latvecs(3,1))))&
          & > epsilon(1.0_dp))) then
        @:RAISE_ERROR(errStatus, -1, "Specified k-value(s) incommensurate with C_n operation.")
      end if
      ! last one is the weight
      ctrl%kWeight(:) = kpts(3, :)
      deallocate(kpts)

    case default
      call detailedError(value1, "Invalid K-point scheme", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine getHelicalKSampling


  !> SCC options that are need for different hamiltonian choices
  subroutine readSccOptions(node, ctrl, geo, errStatus)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    ctrl%tMulliken = .true.

    call getChildValue(node, "ReadInitialCharges", ctrl%tReadChrg, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    if (.not. ctrl%tReadChrg) then
      call getInitialCharges(node, geo, ctrl%initialCharges, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getChildValue(node, "SCCTolerance", ctrl%sccTol, errStatus, 1.0e-5_dp)
    @:PROPAGATE_ERROR(errStatus)

    ! temporararily removed until debugged
    ! call getChildValue(node, "WriteShifts", ctrl%tWriteShifts, errStatus, .false.)
    ! @:PROPAGATE_ERROR(errStatus)
    ctrl%tWriteShifts = .false.

    if (geo%tPeriodic) then
      call getChildValue(node, "EwaldParameter", ctrl%ewaldAlpha, errStatus, 0.0_dp)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(node, "EwaldTolerance", ctrl%tolEwald, errStatus, 1.0e-9_dp)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! self consistency required or not to proceed
    call getChildValue(node, "ConvergentSCCOnly", ctrl%isSccConvRequired, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readSccOptions


  !> Force evaluation options that are need for different hamiltonian choices
  subroutine readForceOptions(node, ctrl, errStatus)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child
    type(string) :: buffer

    call getChildValue(node, "ForceEvaluation", buffer, errStatus, "Traditional", child=child)
    @:PROPAGATE_ERROR(errStatus)
    select case (tolower(unquote(char(buffer))))
    case("traditional")
      ctrl%forceType = forceTypes%orig
    case("dynamicst0")
      ctrl%forceType = forceTypes%dynamicT0
    case("dynamics")
      ctrl%forceType = forceTypes%dynamicTFinite
    case default
      call detailedError(child, "Invalid force evaluation method.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readForceOptions


  !> Options for truncation of the SK data sets at a fixed distance
  subroutine SKTruncations(node, truncationCutOff, skInterMeth, errStatus)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> This is the resulting cutoff distance
    real(dp), intent(out) :: truncationCutOff

    !> Method of the sk interpolation
    integer, intent(in) :: skInterMeth

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    logical :: tHardCutOff
    type(fnode), pointer :: field
    type(string) :: modifier

    ! Artificially truncate the SK table
    call getChildValue(node, "SKMaxDistance", truncationCutOff, errStatus, modifier=modifier,&
        & child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), lengthUnits, field, truncationCutOff, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(node, "HardCutOff", tHardCutOff, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)
    if (tHardCutOff) then
      ! Adjust by the length of the tail appended to the cutoff
      select case(skInterMeth)
      case(skEqGridOld)
        truncationCutOff = truncationCutOff - distFudgeOld
      case(skEqGridNew)
        truncationCutOff = truncationCutOff - distFudge
      end select
    end if
    if (truncationCutOff < epsilon(0.0_dp)) then
      call detailedError(field, "Truncation is shorter than the minimum distance over which SK data&
          & goes to 0", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine SKTruncations


  !> Reads initial charges
  subroutine getInitialCharges(node, geo, initCharges, errStatus)

    !> relevant node in input tree
    type(fnode), pointer :: node

    !> geometry, including atomic type information
    type(TGeometry), intent(in) :: geo

    !> initial atomic charges
    real(dp), allocatable :: initCharges(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child, child2, child3, val
    type(fnodeList), pointer :: children
    integer, allocatable :: pTmpI1(:)
    type(string) :: buffer
    real(dp) :: rTmp
    integer :: ii, jj, iAt

    call getChildValue(node, "InitialCharges", val, errStatus, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true., list=.true.)
    @:PROPAGATE_ERROR(errStatus)

    ! Read either all atom charges, or individual atom specifications
    call getChild(child, "AllAtomCharges", child2, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child2)) then
      allocate(initCharges(geo%nAtom))
      call getChildValue(child2, "", initCharges, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else
      call getChildren(child, "AtomCharge", children)
      if (getLength(children) > 0) then
        allocate(initCharges(geo%nAtom))
        initCharges = 0.0_dp
      end if
      do ii = 1, getLength(children)
        call getItem1(children, ii, child2)
        call getChildValue(child2, "Atoms", buffer, errStatus, child=child3, multiple=.true.)
        @:PROPAGATE_ERROR(errStatus)
        call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species, pTmpI1,&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child2, "ChargePerAtom", rTmp, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        do jj = 1, size(pTmpI1)
          iAt = pTmpI1(jj)
          if (initCharges(iAt) /= 0.0_dp) then
            call detailedWarning(child3, "Previous setting for the charge &
                &of atom" // i2c(iAt) // " overwritten")
          end if
          initCharges(iAt) = rTmp
        end do
        deallocate(pTmpI1)
      end do
      call destroyNodeList(children)
    end if

  end subroutine getInitialCharges


  !> Reads initial spins
  subroutine getInitialSpins(node, geo, nSpin, initSpins, errStatus)

    !> relevant node in input data
    type(fnode), pointer :: node

    !> geometry, including atomic information
    type(TGeometry), intent(in) :: geo

    !> number of spin channels
    integer, intent(in) :: nSpin

    !> initial spins on return
    real(dp), allocatable :: initSpins(:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child, child2, child3, val
    type(fnodeList), pointer :: children
    integer, allocatable :: pTmpI1(:)
    type(string) :: buffer
    real(dp), allocatable :: rTmp(:)
    integer :: ii, jj, iAt

    @:ASSERT(nSpin == 1 .or. nSpin == 3)

    call getChildValue(node, "InitialSpins", val, errStatus, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true., list=.true.)
    @:PROPAGATE_ERROR(errStatus)

    ! Read either all atom spins, or individual spin specifications
    call getChild(child, "AllAtomSpins", child2, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child2)) then
      allocate(initSpins(nSpin, geo%nAtom))
      call getChildValue(child2, "", initSpins, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else
      call getChildren(child, "AtomSpin", children)
      if (getLength(children) > 0) then
        allocate(initSpins(nSpin, geo%nAtom))
        initSpins = 0.0_dp
      end if
      allocate(rTmp(nSpin))
      do ii = 1, getLength(children)
        call getItem1(children, ii, child2)
        call getChildValue(child2, "Atoms", buffer, errStatus, child=child3, multiple=.true.)
        @:PROPAGATE_ERROR(errStatus)
        call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species, pTmpI1,&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child2, "SpinPerAtom", rTmp, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        do jj = 1, size(pTmpI1)
          iAt = pTmpI1(jj)
          if (any(initSpins(:,iAt) /= 0.0_dp)) then
            call detailedWarning(child3, "Previous setting for the spin of atom" // i2c(iAt) //&
                & " overwritten")
          end if
          initSpins(:,iAt) = rTmp
        end do
        deallocate(pTmpI1)
      end do
      deallocate(rTmp)
      call destroyNodeList(children)
    end if

  end subroutine getInitialSpins


  !> Reads numerical differentiation method to be used
  subroutine readDifferentiation(node, ctrl, errStatus)

    !> relevant node in input tree
    type(fnode), pointer, intent(in) :: node

    !> control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> Error status
    type(TStatus), intent(inout) :: errStatus


    !> default of a reasonable choice for round off when using a second order finite difference
    !> formula
    real(dp), parameter :: defDelta = epsilon(1.0_dp)**0.25_dp

    type(string) :: buffer, modifier
    type(fnode), pointer :: val, child

    call getChildValue(node, "Differentiation", val, errStatus, "FiniteDiff", child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(val, buffer)
    select case (char(buffer))
    case ("finitediff")
      ctrl%iDerivMethod = diffTypes%finiteDiff
      call getChildValue(val, "Delta", ctrl%deriv1stDelta, errStatus, defDelta, modifier=modifier,&
          & child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, child, ctrl%deriv1stDelta, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("richardson")
      ctrl%iDerivMethod = diffTypes%richardson
    case default
      call getNodeHSDName(val, buffer)
      call detailedError(child, "Invalid derivative calculation '" // char(buffer) // "'",&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readDifferentiation


  !> Reads the H corrections (H5, Damp)
  subroutine readHCorrection(node, geo, ctrl, errStatus)

    !> Node containing the h-bond correction sub-block.
    type(fnode), pointer, intent(in) :: node

    !> Geometry.
    type(TGeometry), intent(in) :: geo

    !> Control structure
    type(TControl), intent(inout) :: ctrl

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child, child2
    type(string) :: buffer
    real(dp) :: h5ScalingDef
    integer :: iSp

    ! X-H interaction corrections including H5 and damping
    ctrl%tDampH = .false.
    call getChildValue(node, "HCorrection", value1, errStatus, "None", child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)

    select case (char(buffer))

    case ("none")
      ! nothing to do

    case ("damping")
      ! Switch the correction on
      ctrl%tDampH = .true.
      call getChildValue(value1, "Exponent", ctrl%dampExp, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case ("h5")
      allocate(ctrl%h5Input)
      associate (h5Input => ctrl%h5Input)
        call getChildValue(value1, "RScaling", h5Input%rScale, errStatus, 0.714_dp)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(value1, "WScaling", h5Input%wScale, errStatus, 0.25_dp)
        @:PROPAGATE_ERROR(errStatus)
        allocate(h5Input%elementParams(geo%nSpecies))
        call getChild(value1, "H5Scaling", child2, errStatus, requested=.false.)
        @:PROPAGATE_ERROR(errStatus)
        if (.not. associated(child2)) then
          call setChild(value1, "H5scaling", child2, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        do iSp = 1, geo%nSpecies
          select case (geo%speciesNames(iSp))
          case ("O")
            h5ScalingDef = 0.06_dp
          case ("N")
            h5ScalingDef = 0.18_dp
          case ("S")
            h5ScalingDef = 0.21_dp
          case default
            ! Default value is -1, this indicates that the element should be ignored
            h5ScalingDef = -1.0_dp
          end select
          call getChildValue(child2, geo%speciesNames(iSp), h5Input%elementParams(iSp),&
              & errStatus, h5ScalingDef)
          @:PROPAGATE_ERROR(errStatus)
        end do
        h5Input%speciesNames = geo%speciesNames
      end associate

    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Invalid HCorrection '" // char(buffer) // "'", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readHCorrection


  !> Reads Slater-Koster files
  !> Should be replaced with a more sophisticated routine, once the new SK-format has been
  !> established
  subroutine readSKFiles(skFiles, nSpecies, slako, orb, angShells, orbRes, skInterMeth, repPoly,&
      & errStatus, truncationCutOff, rangeSepSK)

    !> List of SK file names to read in for every interaction
    type(TListCharLc), intent(inout) :: skFiles(:,:)

    !> Nr. of species in the system
    integer, intent(in) :: nSpecies

    !> Data type for slako information
    type(TSlater), intent(inout) :: slako

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> For every species, a list of rank one arrays. Each array contains the angular momenta to pick
    !> from the appropriate SK-files.
    type(TListIntR1), intent(inout) :: angShells(:)

    !> Are the Hubbard Us different for each l-shell?
    logical, intent(in) :: orbRes

    !> Method of the sk interpolation
    integer, intent(in) :: skInterMeth

    !> is this a polynomial or spline repulsive?
    logical, intent(in) :: repPoly(:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> Distances to artificially truncate tables of SK integrals
    real(dp), intent(in), optional :: truncationCutOff

    !> if calculation range separated then read omega from end of SK file
    type(TRangeSepSKTag), intent(inout), optional :: rangeSepSK

    integer :: iSp1, iSp2, nSK1, nSK2, iSK1, iSK2, ind, nInteract, iSh1
    integer :: angShell(maxL+1), nShell
    logical :: readRep, readAtomic
    character(lc) :: fileName
    real(dp), allocatable, target :: skHam(:,:), skOver(:,:)
    real(dp) :: dist
    type(TOldSKData), allocatable :: skData12(:,:), skData21(:,:)
    type(TSlakoEqGrid), allocatable :: pSlakoEqGrid1, pSlakoEqGrid2
    type(TSplineRepInp) :: repSplineIn1, repSplineIn2
    type(TPolyRepInp) :: repPolyIn1, repPolyIn2

    ! if artificially cutting the SK tables
    integer :: nEntries

    @:ASSERT(size(skFiles, dim=1) == size(skFiles, dim=2))
    @:ASSERT((size(skFiles, dim=1) > 0) .and. (size(skFiles, dim=1) == nSpecies))
    @:ASSERT(all(shape(repPoly) == shape(skFiles)))

    allocate(slako%skSelf(orb%mShell, nSpecies))
    allocate(slako%skHubbU(orb%mShell, nSpecies))
    allocate(slako%skOcc(orb%mShell, nSpecies))
    allocate(slako%mass(nSpecies))
    slako%skSelf(:,:) = 0.0_dp
    slako%skHubbU(:,:) = 0.0_dp
    slako%skOcc(:,:) = 0.0_dp

    allocate(slako%skHamCont)
    call init(slako%skHamCont, nSpecies)
    allocate(slako%skOverCont)
    call init(slako%skOverCont, nSpecies)
    allocate(slako%pairRepulsives(nSpecies, nSpecies))

    write(stdout, "(A)") "Reading SK-files:"
    lpSp1: do iSp1 = 1, nSpecies
      nSK1 = len(angShells(iSp1))
      lpSp2: do iSp2 = iSp1, nSpecies
        nSK2 = len(angShells(iSp2))
        allocate(skData12(nSK2, nSK1))
        allocate(skData21(nSK1, nSK2))
        ind = 1
        do iSK1 = 1, nSK1
          do iSK2 = 1, nSK2
            readRep = (iSK1 == 1 .and. iSK2 == 1)
            readAtomic = (iSp1 == iSp2 .and. iSK1 == iSK2)
            call get(skFiles(iSp2, iSp1), fileName, ind)
            write(stdOut, "(a)") trim(fileName)
            if (.not. present(rangeSepSK)) then
              if (readRep .and. repPoly(iSp2, iSp1)) then
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, polyRepIn=repPolyIn1)
              elseif (readRep) then
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, iSp1, iSp2,&
                    & splineRepIn=repSplineIn1)
              else
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic)
              end if
            else
              if (readRep .and. repPoly(iSp2, iSp1)) then
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, polyRepIn=repPolyIn1,&
                    & rangeSepSK=rangeSepSK)
              elseif (readRep) then
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, iSp1, iSp2,&
                    & splineRepIn=repSplineIn1, rangeSepSK=rangeSepSK)
              else
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, rangeSepSK=rangeSepSK)
              end if
            end if
            ind = ind + 1
          end do
        end do
        if (iSp1 == iSp2) then
          skData21 = skData12
          if (repPoly(iSp1, iSp2)) then
            repPolyIn2 = repPolyIn1
          else
            repSplineIn2 = repSplineIn1
          end if
          ind = 1
          do iSK1 = 1, nSK1
            call intoArray(angShells(iSp1), angShell, nShell, iSK1)
            do iSh1 = 1, nShell
              slako%skSelf(ind, iSp1) = &
                  &skData12(iSK1,iSK1)%skSelf(angShell(iSh1)+1)
              slako%skOcc(ind:ind, iSp1) = &
                  &skData12(iSK1,iSK1)%skOcc(angShell(iSh1)+1)
              slako%skHubbU(ind, iSp1) = &
                  &skData12(iSK1,iSK1)%skHubbU(angShell(iSh1)+1)
              ind = ind + 1
            end do
          end do
          if (.not. orbRes) then
            slako%skHubbU(2:,iSp1) = slako%skHubbU(1,iSp1)
          end if
          slako%mass(iSp1) = skData12(1,1)%mass
        else
          ind = 1
          do iSK2 = 1, nSK2
            do iSK1 = 1, nSK1
              readRep = (iSK1 == 1 .and. iSK2 == 1)
              call get(sKFiles(iSp1, iSp2), fileName, ind)
              if (readRep .and. repPoly(iSp1, iSp2)) then
                call readFromFile(skData21(iSK1,iSK2), fileName, readAtomic, &
                    &polyRepIn=repPolyIn2)
              elseif (readRep) then
                call readFromFile(skData21(iSK1,iSK2), fileName, readAtomic, &
                    &iSp2, iSp1, splineRepIn=repSplineIn2)
              else
                call readFromFile(skData21(iSK1,iSK2), fileName, readAtomic)
              end if
              ind = ind + 1
            end do
          end do
        end if

        ! Check for SK and repulsive consistentcy
        call checkSKCompElec(skData12, skData21, iSp1, iSp2, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        if (repPoly(iSp1, iSp2)) then
          call checkSKCompRepPoly(repPolyIn1, repPolyIn2, iSp1, iSp2, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        else
          call checkSKCompRepSpline(repSplineIn1, repSplineIn2, iSp1, iSp2, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if

        ! Create full H/S table for all interactions of iSp1-iSp2
        nInteract = getNSKIntegrals(iSp1, iSp2, orb)
        allocate(skHam(size(skData12(1,1)%skHam, dim=1), nInteract))
        allocate(skOver(size(skData12(1,1)%skOver, dim=1), nInteract))
        call getFullTable(skHam, skOver, skData12, skData21, angShells(iSp1),&
            & angShells(iSp2))

        ! Add H/S tables to the containers for iSp1-iSp2
        dist = skData12(1,1)%dist
        if (present(truncationCutOff)) then
          nEntries = floor(truncationCutOff / dist)
          nEntries = min(nEntries, size(skData12(1,1)%skHam, dim=1))
        else
          nEntries = size(skData12(1,1)%skHam, dim=1)
        end if
        allocate(pSlakoEqGrid1, pSlakoEqGrid2)
        call init(pSlakoEqGrid1, dist, skHam(:nEntries,:), skInterMeth)
        call init(pSlakoEqGrid2, dist, skOver(:nEntries,:), skInterMeth)
        call addTable(slako%skHamCont, pSlakoEqGrid1, iSp1, iSp2)
        call addTable(slako%skOverCont, pSlakoEqGrid2, iSp1, iSp2)
        deallocate(skHam)
        deallocate(skOver)
        if (iSp1 /= iSp2) then
          ! Heteronuclear interactions: the same for the reverse interaction
          allocate(skHam(size(skData12(1,1)%skHam, dim=1), nInteract))
          allocate(skOver(size(skData12(1,1)%skOver, dim=1), nInteract))
          call getFullTable(skHam, skOver, skData21, skData12, angShells(iSp2),&
              & angShells(iSp1))
          allocate(pSlakoEqGrid1, pSlakoEqGrid2)
          call init(pSlakoEqGrid1, dist, skHam(:nEntries,:), skInterMeth)
          call init(pSlakoEqGrid2, dist, skOver(:nEntries,:), skInterMeth)
          call addTable(slako%skHamCont, pSlakoEqGrid1, iSp2, iSp1)
          call addTable(slako%skOverCont, pSlakoEqGrid2, iSp2, iSp1)
          deallocate(skHam)
          deallocate(skOver)
        end if
        deallocate(skData12)
        deallocate(skData21)

        ! Create repulsive container

        ! Add repulsives to the containers.
        if (repPoly(iSp2, iSp1)) then
          slako%pairRepulsives(iSp2, iSp1)%item = TPolyRep(repPolyIn1)
        else
          slako%pairRepulsives(iSp2, iSp1)%item = TSplineRep(repSplineIn1)
          deallocate(repSplineIn1%xStart)
          deallocate(repSplineIn1%spCoeffs)
        end if
        if (iSp1 /= iSp2) then
          if (repPoly(iSp1, iSp2)) then
            slako%pairRepulsives(iSp1, iSp2)%item = TPolyRep(repPolyIn2)
          else
            slako%pairRepulsives(iSp1, iSp2)%item = TSplineRep(repSplineIn2)
            deallocate(repSplineIn2%xStart)
            deallocate(repSplineIn2%spCoeffs)
          end if
        end if
      end do lpSp2
    end do lpSp1
    write(stdout, "(A)") "Done."

  end subroutine readSKFiles


  !> Checks if the provided set of SK-tables for a the interactions A-B and B-A are consistent
  subroutine checkSKCompElec(skData12, skData21, sp1, sp2, errStatus)

    !> Slater-Koster integral set for the interaction A-B
    type(TOldSKData), intent(in), target :: skData12(:,:)

    !> Slater-Koster integral set for the interaction B-A
    type(TOldSKData), intent(in), target :: skData21(:,:)

    !> Species number for A (for error messages)
    integer, intent(in) :: sp1

    !> Species number for B (for error messages)
    integer, intent(in) :: sp2

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: iSK1, iSK2, nSK1, nSK2
    integer :: nGrid
    real(dp) :: dist
    type(TOldSKData), pointer :: pSK12, pSK21
    character(lc) :: errorStr

    nSK1 = size(skData12, dim=2)
    nSK2 = size(skData12, dim=1)

    @:ASSERT(size(skData21, dim=1) == nSK1)
    @:ASSERT(size(skData21, dim=2) == nSK2)

    nGrid = skData12(1,1)%nGrid
    dist = skData12(1,1)%dist

    ! All SK files should have the same grid separation and table length
    nGrid = skData12(1,1)%nGrid
    dist = skData12(1,1)%dist
    do iSK1 = 1, nSK1
      do iSK2 = 1, nSK2
        pSK12 => skData12(iSK2, iSK1)
        pSK21 => skData21(iSK1, iSK2)

        if (pSK12%dist /= dist .or. pSK21%dist /= dist) then
          write (errorStr, "(A,I2,A,I2)") "Incompatible SK grid separations &
              &for species ", sp1, ", ", sp2
          @:RAISE_ERROR(errStatus, -1, errorStr)
        end if
        if (pSK12%nGrid /= nGrid .or. pSK21%nGrid /= nGrid) then
          write (errorStr, "(A,I2,A,I2)") "Incompatible SK grid lengths for &
              &species pair ", sp1, ", ", sp2
          @:RAISE_ERROR(errStatus, -1, errorStr)
        end if
      end do
    end do

  end subroutine checkSKCompElec


  !> Checks if the provided repulsive splines for A-B and B-A are compatible
  subroutine checkSKCompRepSpline(repIn1, repIn2, sp1, sp2, errStatus)

    !> Repulsive spline for interaction A-B
    type(TSplineRepInp), intent(in) :: repIn1

    !> Repulsive spline for interaction B-A
    type(TSplineRepInp), intent(in) :: repIn2

    !> Number of species A (for error messages only)
    integer, intent(in) :: sp1

    !> Number of species B (for error messages only)
    integer, intent(in) :: sp2

    !> Error status
    type(TStatus), intent(inout) :: errStatus


    !> Tolerance for the agreement in the repulsive data
    real(dp), parameter :: tolRep = 1.0e-8_dp


    !> string for error return
    character(lc) :: errorStr

    ! Repulsives for A-B and B-A should be the same
    if (size(repIn1%xStart) /= size(repIn2%xStart)) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible nr. of repulsive &
          &intervals for species pair ", sp1, "-", sp2, "."
      @:RAISE_ERROR(errStatus, -1, errorStr)
    end if
    if (maxval(abs(repIn1%xStart - repIn2%xStart)) > tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &intervals for species pair ", sp1, "-", sp2, "."
      @:RAISE_ERROR(errStatus, -1, errorStr)
    end if
    if (maxval(abs(repIn1%spCoeffs - repIn2%spCoeffs)) > tolRep &
        &.or. maxval(abs(repIn1%spLastCoeffs - repIn2%spLastCoeffs)) &
        &> tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &coefficients for species pair ", sp1, "-", sp2, "."
      @:RAISE_ERROR(errStatus, -1, errorStr)
    end if
    if (maxval(abs(repIn1%expCoeffs - repIn2%expCoeffs)) > tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &exp. coefficients for species pair ", sp1, "-", sp2, "."
      @:RAISE_ERROR(errStatus, -1, errorStr)
    end if
    if (abs(repIn1%cutoff - repIn2%cutoff) > tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &cutoffs for species pair ", sp1, "-", sp2, "."
      @:RAISE_ERROR(errStatus, -1, errorStr)
    end if

  end subroutine checkSKCompRepSpline


  !> Checks if repulsive polynomials for A-B and B-A are compatible
  subroutine checkSKCompRepPoly(repIn1, repIn2, sp1, sp2, errStatus)

    !> Repulsive polynomial for interaction A-B
    type(TPolyRepInp), intent(in) :: repIn1

    !> Repulsive polynomial for interaction B-A
    type(TPolyRepInp), intent(in) :: repIn2

    !> Number of species A (for error messages only)
    integer, intent(in) :: sp1

    !> Number of species B (for error messages only)
    integer, intent(in) :: sp2

    !> Error status
    type(TStatus), intent(inout) :: errStatus


    !> for error string return
    character(lc) :: errorStr

    if (any(repIn1%polyCoeffs /= repIn2%polyCoeffs)) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive polynomial &
          &coefficients  for the species pair ", sp1, "-", sp2, "."
      @:RAISE_ERROR(errStatus, -1, errorStr)
    end if
    if (repIn1%cutoff /= repIn2%cutoff) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive cutoffs  &
          &for the species pair ", sp1, "-", sp2, "."
      @:RAISE_ERROR(errStatus, -1, errorStr)
    end if

  end subroutine checkSKCompRepPoly


  !> Returns the nr. of Slater-Koster integrals necessary to describe the interactions between two
  !> species
  pure function getNSKIntegrals(sp1, sp2, orb) result(nInteract)

    !> Index of the first species
    integer, intent(in) :: sp1

    !> Index of the second species
    integer, intent(in) :: sp2

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Nr. of Slater-Koster interactions
    integer :: nInteract

    integer :: iSh1, iSh2

    nInteract = 0
    do iSh1 = 1, orb%nShell(sp1)
      do iSh2 = 1, orb%nShell(sp2)
        nInteract = nInteract + min(orb%angShell(iSh2, sp2), orb%angShell(iSh1, sp1)) + 1
      end do
    end do

  end function getNSKIntegrals


  !> Creates from the columns of the Slater-Koster files for A-B and B-A a full table for A-B,
  !> containing all integrals.
  subroutine getFullTable(skHam, skOver, skData12, skData21, angShells1, angShells2)

    !> Resulting table of H integrals
    real(dp), intent(out) :: skHam(:,:)

    !> Resulting table of S integrals
    real(dp), intent(out) :: skOver(:,:)

    !> Contains all SK files describing interactions for A-B
    type(TOldSKData), intent(in), target :: skData12(:,:)

    !> Contains all SK files describing interactions for B-A
    type(TOldSKData), intent(in), target :: skData21(:,:)

    !> Angular momenta to pick from the SK-files for species A
    type(TListIntR1), intent(inout) :: angShells1

    !> Angular momenta to pick from the SK-files for species B
    type(TListIntR1), intent(inout) :: angShells2

    integer :: ind, iSK1, iSK2, iSh1, iSh2, nSh1, nSh2, l1, l2, lMin, lMax, mm
    integer :: angShell1(maxL+1), angShell2(maxL+1)
    real(dp), pointer :: pHam(:,:), pOver(:,:)


    !> Maps (mm, l1, l2 ) onto an element in the SK table.
    !> l2 >= l1 (l1 = 0, 1, ...; l2 = 0, 1, ...), m <= l1.
    integer, parameter :: skMap(0:maxL, 0:maxL, 0:maxL) &
        &= reshape((/&
        &20, 0,  0,  0,  19,  0,  0,  0,  18,  0,  0,  0,  17,  0,  0,  0,&
        & 0, 0,  0,  0,  15, 16,  0,  0,  13, 14,  0,  0,  11, 12,  0,  0,&
        & 0, 0,  0,  0,   0,  0,  0,  0,   8,  9, 10,  0,   5,  6,  7,  0,&
        & 0, 0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   1,  2,  3,  4/),&
        &(/maxL + 1, maxL + 1, maxL + 1/))

    ind = 1
    do iSK1 = 1, len(angShells1)
      call intoArray(angShells1, angShell1, nSh1, iSK1)
      do iSh1 = 1, nSh1
        l1 = angShell1(iSh1)
        do iSK2 = 1, len(angShells2)
          call intoArray(angShells2, angShell2, nSh2, iSK2)
          do iSh2 = 1, nSh2
            l2 = angShell2(iSh2)
            if (l1 <= l2) then
              pHam => skData12(iSK2,iSK1)%skHam
              pOver => skData12(iSK2,iSK1)%skOver
              lMin = l1
              lMax = l2
            else
              pHam => skData21(iSK1,iSK2)%skHam
              pOver => skData21(iSK1,iSK2)%skOver
              lMin = l2
              lMax = l1
            end if
            do mm = 0, lMin
              ! Safety check, if array size are appropriate
              @:ASSERT(all(shape(skHam) >= (/ size(pHam, dim=1), ind /)))
              @:ASSERT(all(shape(skOver) >= (/ size(pOver, dim=1), ind /)))
              @:ASSERT(size(pHam, dim=1) == size(pOver, dim=1))
              skHam(:,ind) = pHam(:,skMap(mm,lMax,lMin))
              skOver(:,ind) = pOver(:,skMap(mm,lMax,lMin))
              ind = ind + 1
            end do
          end do
        end do
      end do
    end do

  end subroutine getFullTable


  !> Reads the option block
  subroutine readOptions(node, ctrl, errStatus)

    !> Node to parse
    type(fnode), pointer :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child
    type(string) :: strBuffer
    logical :: tWriteDetailedOutDef

  #:if WITH_SOCKETS
    tWriteDetailedOutDef = .not. allocated(ctrl%socketInput)
  #:else
    tWriteDetailedOutDef = .true.
  #:endif
    call getChildValue(node, "WriteDetailedOut", ctrl%tWriteDetailedOut, errStatus,&
        & tWriteDetailedOutDef)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(node, "WriteAutotestTag", ctrl%tWriteTagged, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "WriteDetailedXML", ctrl%tWriteDetailedXML, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "WriteResultsTag", ctrl%tWriteResultsTag, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)

    if (.not.(ctrl%tMD.or.ctrl%isGeoOpt.or.allocated(ctrl%geoOpt))) then
      if (ctrl%tSCC) then
        call getChildValue(node, "RestartFrequency", ctrl%restartFreq, errStatus, 20)
        @:PROPAGATE_ERROR(errStatus)
      else
        ctrl%restartFreq = 0
      end if
    end if
    call getChildValue(node, "RandomSeed", ctrl%iSeed, errStatus, 0, child=child)
    @:PROPAGATE_ERROR(errStatus)
    if (ctrl%iSeed < 0) then
      call detailedError(child, "Random seed must be greater or equal zero", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    call getChildValue(node, "WriteHS", ctrl%tWriteHS, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "WriteRealHS", ctrl%tWriteRealHS, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call renameChildren(node, "MinimizeMemoryUsage", "MinimiseMemoryUsage")
    call getChildValue(node, "MinimiseMemoryUsage", ctrl%tMinMemory, errStatus, .false.,&
        & child=child)
    @:PROPAGATE_ERROR(errStatus)
    if (ctrl%tMinMemory) then
      call detailedWarning(child, "Memory minimisation is not working currently, normal calculation&
          & will be used instead")
    end if
    call getChildValue(node, "ShowFoldedCoords", ctrl%tShowFoldedCoord, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
  #:if DEBUG > 0
    call getChildValue(node, "TimingVerbosity", ctrl%timingLevel, errStatus, -1)
    @:PROPAGATE_ERROR(errStatus)
  #:else
    call getChildValue(node, "TimingVerbosity", ctrl%timingLevel, errStatus, 1)
    @:PROPAGATE_ERROR(errStatus)
  #:endif

    if (ctrl%tReadChrg) then
      call getChildValue(node, "ReadChargesAsText", ctrl%tReadChrgAscii, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getChildValue(node, "WriteCharges", ctrl%tWriteCharges, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)
    if (ctrl%tWriteCharges) then
      call getChildValue(node, "WriteChargesAsText", ctrl%tWriteChrgAscii, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ctrl%tSkipChrgChecksum = .false.
    if (.not. ctrl%tFixEf .and. ctrl%tReadChrg) then
      call getChildValue(node, "SkipChargeTest", ctrl%tSkipChrgChecksum, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call readBinaryAccessTypes(node, ctrl%binaryAccessTypes, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readOptions


  !> Reads in dispersion related settings
  subroutine readDispersion(node, geo, input, nrChrg, tSCC, errStatus)

    !> Node to parse
    type(fnode), pointer :: node

    !> geometry, including atomic information
    type(TGeometry), intent(in) :: geo

    !> dispersion data on exit
    type(TDispersionInp), intent(out) :: input

    !> net charge
    real(dp), intent(in) :: nrChrg

    !> SCC calculation?
    logical, intent(in) :: tScc

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: dispModel
    type(string) :: buffer

    call getChildValue(node, "", dispModel, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(dispModel, buffer)
    select case (char(buffer))
    case ("slaterkirkwood")
      allocate(input%slakirk)
      call readDispSlaKirk(dispModel, geo, input%slakirk, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("lennardjones")
      allocate(input%uff)
      call readDispVdWUFF(dispModel, geo, input%uff, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("dftd3")
      allocate(input%dftd3)
      call readDFTD3(dispModel, geo, input%dftd3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("simpledftd3")
      allocate(input%sdftd3)
      call readSimpleDFTD3(dispModel, geo, input%sdftd3, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("dftd4")
      allocate(input%dftd4)
      call readDispDFTD4(dispModel, geo, input%dftd4, nrChrg, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("ts")
  #:if WITH_MBD
      allocate(input%mbd)
      call readDispTs(dispModel, input%mbd, errStatus)
      @:PROPAGATE_ERROR(errStatus)
  #:else
      call detailedError(node, "Program must be compiled with the mbd library for TS-dispersion",&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
  #:endif
    case ("mbd")
  #:if WITH_MBD
      allocate(input%mbd)
      call readDispMbd(dispModel, input%mbd, errStatus)
      @:PROPAGATE_ERROR(errStatus)
  #:else
      call detailedError(node, "Program must be compiled with the mbd library for MBD-dispersion",&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
  #:endif
    case default
      call detailedError(node, "Invalid dispersion model name.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readDispersion


  !> Reads in the dispersion input data for the Slater-Kirkwood dispersion model
  subroutine readDispSlaKirk(node, geo, input, errStatus)

    !> Node to process
    type(fnode), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the dispersion module on exit
    type(TDispSlaKirkInp), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, value2, child, child2, child3
    type(string) :: buffer, modifier, modifier2, modifiers(3)
    real(dp), allocatable :: tmpR2(:,:), tmp2R2(:,:), rCutoffs(:)
    real(dp) :: mCutoff, rTmp
    integer :: iAt1, iAt2f, iSp1, iSp2, iNeigh
    integer, allocatable :: nNeighs(:)
    real(dp), allocatable :: cellVec(:,:), rCellVec(:,:)
    real(dp), allocatable :: coords(:,:)
    integer, allocatable :: img2CentCell(:), iCellVec(:)
    integer :: nAllAtom
    type(TNeighbourList) :: neighs

    allocate(tmpR2(3, geo%nAtom))
    allocate(input%polar(geo%nAtom))
    allocate(input%rWaals(geo%nAtom))
    allocate(input%charges(geo%nAtom))
    call getChildValue(node, "PolarRadiusCharge", value1, errStatus, child=child, modifier=modifier)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    select case (char(buffer))
    case (textNodeName)
      call getChildValue(child, "", tmpR2, errStatus, modifier=modifier)
      @:PROPAGATE_ERROR(errStatus)
      if (len(modifier) > 0) then
        call splitModifier(char(modifier), child, modifiers, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(modifiers(1)), volumeUnits, child, tmpR2(1,:),&
            & errStatus, .false.)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(modifiers(2)), lengthUnits, child, tmpR2(2,:),&
            & errStatus, .false.)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(modifiers(3)), chargeUnits, child, tmpR2(3,:),&
            & errStatus, .false.)
        @:PROPAGATE_ERROR(errStatus)
      end if

    case ("hybriddependentpol")
      if (len(modifier) > 0) then
        call detailedError(child, "PolarRadiusCharge is not allowed to carry &
            &a modifier, if the HybridDependentPol method is used.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      allocate(rCutoffs(geo%nSpecies))
      allocate(tmp2R2(13, geo%nSpecies))
      do iSp1 = 1, geo%nSpecies
        call getChildValue(value1, geo%speciesNames(iSp1), value2, errStatus, child=child2,&
            & dummyValue=.true.)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child2, "CovalentRadius", rCutoffs(iSp1), errStatus, modifier=modifier2,&
            & child=child3)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(modifier2), lengthUnits, child3, rCutoffs(iSp1), errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call renameChildren(child2, "HybridPolarizations", "HybridPolarisations")
        call getChildValue(child2, "HybridPolarisations", tmp2R2(:, iSp1), errStatus,&
            & modifier=modifier2, child=child3)
        @:PROPAGATE_ERROR(errStatus)
        if (len(modifier2) > 0) then
          call splitModifier(char(modifier2), child, modifiers, errStatus)
          @:PROPAGATE_ERROR(errStatus)
          call convertUnitHsd(char(modifiers(1)), volumeUnits, child, tmp2R2(1:6, iSp1), errStatus,&
              & .false.)
          @:PROPAGATE_ERROR(errStatus)
          call convertUnitHsd(char(modifiers(2)), lengthUnits, child, tmp2R2(7:12, iSp1),&
              & errStatus, .false.)
          @:PROPAGATE_ERROR(errStatus)
          call convertUnitHsd(char(modifiers(3)), chargeUnits, child, tmp2R2(13, iSp1), errStatus,&
              & .false.)
          @:PROPAGATE_ERROR(errStatus)
        end if
      end do
      mCutoff = 2.0_dp * maxval(rCutoffs)
      if (geo%tPeriodic) then
        call getCellTranslations(cellVec, rCellVec, geo%latVecs, geo%recVecs2p, mCutoff)
      else
        allocate(cellVec(3, 1))
        allocate(rCellVec(3, 1))
        cellVec(:, 1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
        rCellVec(:, 1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
      end if
      call TNeighbourlist_init(neighs, geo%nAtom, 10)
      if (geo%tPeriodic) then
        ! Make some guess for the nr. of all interacting atoms
        nAllAtom = int((real(geo%nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
      else
        nAllAtom = geo%nAtom
      end if
      allocate(coords(3, nAllAtom))
      allocate(img2CentCell(nAllAtom))
      allocate(iCellVec(nAllAtom))
      call updateNeighbourList(coords, img2CentCell, iCellVec, neighs, nAllAtom, geo%coords,&
          & mCutoff, rCellVec, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      allocate(nNeighs(geo%nAtom))
      nNeighs(:) = 0
      do iAt1 = 1, geo%nAtom
        iSp1 = geo%species(iAt1)
        do iNeigh = 1, neighs%nNeighbour(iAt1)
          iAt2f = img2CentCell(neighs%iNeighbour(iNeigh, iAt1))
          iSp2 = geo%species(iAt2f)
          rTmp = rCutoffs(iSp1) + rCutoffs(iSp2)
          if (neighs%neighDist2(iNeigh, iAt1) <= rTmp**2) then
            nNeighs(iAt1) = nNeighs(iAt1) + 1
            nNeighs(iAt2f) = nNeighs(iAt2f) + 1
          end if
        end do
      end do
      do iAt1 = 1, geo%nAtom
        iSp1 = geo%species(iAt1)
        if (nNeighs(iAt1) <= 4 ) then
          tmpR2(1, iAt1) = tmp2R2(1+nNeighs(iAt1), iSp1)
          tmpR2(2, iAt1) = tmp2R2(7+nNeighs(iAt1), iSp1)
        else
          tmpR2(1, iAt1) = tmp2R2(6, iSp1)
          tmpR2(2, iAt1) = tmp2R2(12, iSp1)
        end if
        tmpR2(3, iAt1) = tmp2R2(13, iSp1)
      end do

    case default
      call detailedError(value1, "Invalid method for PolarRadiusCharge.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    input%polar(:) = tmpR2(1,:)
    input%rWaals(:) = tmpR2(2,:)
    input%charges(:) = tmpR2(3,:)

  end subroutine readDispSlaKirk


  !> Reads in initialization data for the UFF dispersion model
  subroutine readDispVdWUFF(node, geo, input, errStatus)

    !> Node to process
    type(fnode), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit
    type(TDispUffInp), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(string) :: buffer
    type(fnode), pointer :: child, value1, child2
    integer :: iSp
    logical :: found

    call getChildValue(node, "Parameters", value1, errStatus, child=child)
    @:PROPAGATE_ERROR(errStatus)
    allocate(input%distances(geo%nSpecies))
    allocate(input%energies(geo%nSpecies))
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case("uffparameters")
      do iSp = 1, geo%nSpecies
        call getUffValues(geo%speciesNames(iSp), input%distances(iSp), &
            &input%energies(iSp), found)
        if (.not. found) then
          call detailedError(value1, "UFF parameters for species '" // geo%speciesNames(iSp)&
              & // "' not found.", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      end do
    case default
      call setUnprocessed(value1)
      do iSp = 1, geo%nSpecies
        call getChild(child, geo%speciesNames(iSp), child2, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child2, "Distance", input%distances(iSp), errStatus, modifier=buffer)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(buffer), lengthUnits, child, input%distances(iSp), errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child2, "Energy", input%energies(iSp), errStatus, modifier=buffer)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(buffer), energyUnits, child, input%energies(iSp), errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end do
    end select

  end subroutine readDispVdWUFF


  !> Reads in initialization data for the DFTD3 dispersion module
  subroutine readDFTD3(node, geo, input, errStatus)

    !> Node to process.
    type(fnode), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit.
    type(TSDFTD3Input), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: iSp
    integer, allocatable :: izpDefault(:)
    type(fnode), pointer :: child, childval
    type(string) :: buffer
    integer, parameter :: d3MaxNum = 94
    logical :: unknownSpecies, threebody

    call getChildValue(node, "Damping", childval, errStatus, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(childval, buffer)
    select case (char(buffer))
    case ("beckejohnson")
      input%dampingFunction = dampingFunction%rational
      call getChildValue(childval, "a1", input%a1, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(childval, "a2", input%a2, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("zerodamping")
      input%dampingFunction = dampingFunction%zero
      call getChildValue(childval, "sr6", input%sr6, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(childval, "alpha6", input%alpha6, errStatus, default=14.0_dp)
      @:PROPAGATE_ERROR(errStatus)
    case ("modifiedzerodamping")
      input%dampingFunction = dampingFunction%mzero
      call getChildValue(childval, "sr6", input%sr6, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(childval, "beta", input%beta, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(childval, "alpha6", input%alpha6, errStatus, default=14.0_dp)
      @:PROPAGATE_ERROR(errStatus)
    case default
      call getNodeHSDName(childval, buffer)
      call detailedError(child, "Invalid damping method '" // char(buffer) // "'", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select
    call getChildValue(node, "s6", input%s6, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "s8", input%s8, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "cutoff", input%cutoff, errStatus, default=sqrt(9000.0_dp),&
        & modifier=buffer, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoff, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "cutoffcn", input%cutoffCN, errStatus, default=40.0_dp,&
        & modifier=buffer, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoffCN, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "threebody", threebody, errStatus, default=.false.)
    @:PROPAGATE_ERROR(errStatus)
    input%s9 = merge(1.0_dp, 0.0_dp, threebody)
    ! D3H5 - additional H-H repulsion
    call getChildValue(node, "hhrepulsion", input%hhrepulsion, errStatus, default=.false.)
    @:PROPAGATE_ERROR(errStatus)

    ! Initialize default atomic numbers
    allocate(izpDefault(size(geo%speciesNames)))
    do iSp = 1, size(geo%speciesNames)
      izpDefault(iSp) = symbolToNumber(geo%speciesNames(iSp))
    end do

    ! See if we find user specified overwrites for atomic numbers
    call getChild(node, "AtomicNumbers", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      allocate(input%izp(size(geo%speciesNames)))
      call readSpeciesList(child, geo%speciesNames, input%izp, errStatus, izpDefault)
      @:PROPAGATE_ERROR(errStatus)
      deallocate(izpDefault)
    else
      call move_alloc(izpDefault, input%izp)
    end if

    unknownSpecies = .false.
    do iSp = 1, size(geo%speciesNames)
      if (input%izp(iSp) <= 0 .or. input%izp(iSp) > d3MaxNum) then
        unknownSpecies = .true.
        call warning("Species '"//trim(geo%speciesNames(iSp))// &
          & "' is not supported by DFT-D3")
      end if
    end do
    if (unknownSpecies) then
      call detailedError(node, "DFT-D3 does not support all species present", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine readDFTD3


  !> Reads in initialization data for the simple D3 dispersion model.
  subroutine readSimpleDFTD3(node, geo, input, errStatus)

    !> Node to process.
    type(fnode), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit.
    type(TSimpleDftD3Input), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child
    type(string) :: buffer

    call getChildValue(node, "s6", input%s6, errStatus, default=1.0_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "s8", input%s8, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "s10", input%s10, errStatus, default=0.0_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "a1", input%a1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "a2", input%a2, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "alpha", input%alpha, errStatus, default=14.0_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "weightingFactor", input%weightingFactor, errStatus, default=4.0_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "cutoffInter", input%cutoffInter, errStatus, default=64.0_dp,&
        & modifier=buffer, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoffInter, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call readCoordinationNumber(node, input%cnInput, geo, "exp", 0.0_dp, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readSimpleDFTD3


  !> Reads in initialization data for the D4 dispersion model.
  !>
  !> The D4 dispersion model is usually constructed in a failsafe way, so
  !> it only requires to know the damping parameters s8, a1 and a2.
  !> Here we additionally require a s9, since the non-addititive contributions
  !> tend to be expensive especially in the tight-binding context, s9 = 0.0_dp
  !> will disable the calculation.
  subroutine readDispDFTD4(node, geo, input, nrChrg, errStatus)

    !> Node to process.
    type(fnode), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit.
    type(TDispDftD4Inp), intent(out) :: input

    !> Net charge of the system.
    real(dp), intent(in) :: nrChrg

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: iSp
    integer, allocatable :: izpDefault(:)
    type(fnode), pointer :: value1, child
    type(string) :: buffer
    real(dp), allocatable :: d4Chi(:), d4Gam(:), d4Kcn(:), d4Rad(:)
    integer, parameter :: d4MaxNum = 86
    logical :: unknownSpecies

    call getChildValue(node, "s6", input%s6, errStatus, default=1.0_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "s8", input%s8, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "s9", input%s9, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "s10", input%s10, errStatus, default=0.0_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "a1", input%a1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "a2", input%a2, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "alpha", input%alpha, errStatus, default=16.0_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "WeightingFactor", input%weightingFactor, errStatus, default=6.0_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "ChargeSteepness", input%chargeSteepness, errStatus, default=2.0_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "ChargeScale", input%chargeScale, errStatus, default=3.0_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "CutoffInter", input%cutoffInter, errStatus, default=64.0_dp,&
        & modifier=buffer, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoffInter, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "CutoffThree", input%cutoffThree, errStatus, default=40.0_dp,&
        & modifier=buffer, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoffThree, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(node, "ChargeModel", value1, errStatus, "EEQ", child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(value1, "Unknown method '"//char(buffer)//"' for ChargeModel", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("selfconsistent")
      input%selfConsistent = .true.
    case ("eeq")
      allocate(input%eeqInput)
      allocate(d4Chi(geo%nSpecies))
      d4Chi(:) = getEeqChi(geo%speciesNames)
      allocate(d4Gam(geo%nSpecies))
      d4Gam(:) = getEeqGam(geo%speciesNames)
      allocate(d4Kcn(geo%nSpecies))
      d4Kcn(:) = getEeqKcn(geo%speciesNames)
      allocate(d4Rad(geo%nSpecies))
      d4Rad(:) = getEeqRad(geo%speciesNames)
      call readEeqModel(value1, input%eeqInput, geo, nrChrg, d4Chi, d4Gam, d4Kcn, d4Rad, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    ! Initialize default atomic numbers
    allocate(izpDefault(size(geo%speciesNames)))
    do iSp = 1, size(geo%speciesNames)
      izpDefault(iSp) = symbolToNumber(geo%speciesNames(iSp))
    end do

    ! See if we find user specified overwrites for atomic numbers
    call getChild(node, "AtomicNumbers", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      allocate(input%izp(size(geo%speciesNames)))
      call readSpeciesList(child, geo%speciesNames, input%izp, errStatus, izpDefault)
      @:PROPAGATE_ERROR(errStatus)
      deallocate(izpDefault)
    else
      call move_alloc(izpDefault, input%izp)
    end if

    call readCoordinationNumber(node, input%cnInput, geo, "Cov", 0.0_dp, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    unknownSpecies = .false.
    do iSp = 1, size(geo%speciesNames)
      if (input%izp(iSp) <= 0 .or. input%izp(iSp) > d4MaxNum) then
        unknownSpecies = .true.
        call warning("Species '"//trim(geo%speciesNames(iSp))// &
          & "' is not supported by DFT-D4")
      end if
    end do
    if (unknownSpecies) then
      call detailedError(node, "DFT-D4 does not support all species present", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine readDispDFTD4


  !> Read settings regarding the EEQ charge model
  subroutine readEeqModel(node, input, geo, nrChrg, kChiDefault, kGamDefault, kKcnDefault,&
      & kRadDefault, errStatus)

    !> Node to process.
    type(fnode), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit.
    type(TEeqInput), intent(out) :: input

    !> Net charge of the system.
    real(dp), intent(in) :: nrChrg

    !> Electronegativities default values
    real(dp), intent(in) :: kChiDefault(:)

    !> Chemical hardnesses default values
    real(dp), intent(in) :: kGamDefault(:)

    !> CN scaling default values
    real(dp), intent(in) :: kKcnDefault(:)

    !> Charge widths default values
    real(dp), intent(in) :: kRadDefault(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child
    type(string) :: buffer

    input%nrChrg = nrChrg

    allocate(input%chi(geo%nSpecies))
    allocate(input%gam(geo%nSpecies))
    allocate(input%kcn(geo%nSpecies))
    allocate(input%rad(geo%nSpecies))

    call getChildValue(node, "Chi", value1, errStatus, "Defaults", child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' for chi", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%chi, errStatus, kChiDefault)
      @:PROPAGATE_ERROR(errStatus)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%chi, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    call getChildValue(node, "Gam", value1, errStatus, "Defaults", child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' for gam", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%gam, errStatus, kGamDefault)
      @:PROPAGATE_ERROR(errStatus)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%gam, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    call getChildValue(node, "Kcn", value1, errStatus, "Defaults", child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' for kcn", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%kcn, errStatus, kKcnDefault)
      @:PROPAGATE_ERROR(errStatus)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%kcn, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    call getChildValue(node, "Rad", value1, errStatus, "Defaults", child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' for rad", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%rad, errStatus, kRadDefault)
      @:PROPAGATE_ERROR(errStatus)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%rad, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    call getChildValue(node, "Cutoff", input%cutoff, errStatus, default=40.0_dp, modifier=buffer,&
        & child=child)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoff, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(node, "EwaldParameter", input%parEwald, errStatus, 0.0_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "EwaldTolerance", input%tolEwald, errStatus, 1.0e-9_dp)
    @:PROPAGATE_ERROR(errStatus)

    call readCoordinationNumber(node, input%cnInput, geo, "Erf", 8.0_dp, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readEeqModel


  !> Read in coordination number settings
  subroutine readCoordinationNumber(node, input, geo, cnDefault, cutDefault, errStatus)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TCNInput), intent(inout) :: input

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Default value for the coordination number type
    character(len=*), intent(in) :: cnDefault

    !> Default value for the maximum CN used for cutting (0 turns it off)
    real(dp), intent(in) :: cutDefault

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, value2, child, child2, field
    type(string) :: buffer, modifier
    real(dp), allocatable :: kENDefault(:), kRadDefault(:)

    call getChildValue(node, "CoordinationNumber", value1, errStatus, cnDefault, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)

    select case(char(buffer))
    case default
      call detailedError(child, "Invalid coordination number type specified", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case("exp")
      input%cnType = cnType%exp
    case("erf")
      input%cnType = cnType%erf
    case("cov")
      input%cnType = cnType%cov
    case("gfn")
      input%cnType = cnType%gfn
    end select

    call getChildValue(value1, "CutCN", input%maxCN, errStatus, cutDefault, child=child2)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(value1, "Cutoff", input%rCutoff, errStatus, 40.0_dp, modifier=modifier,&
        & child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), lengthUnits, field, input%rCutoff, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    allocate(input%en(geo%nSpecies))
    if (input%cnType == cnType%cov) then
      call getChildValue(value1, "Electronegativities", value2, errStatus, "PaulingEN",&
          & child=child2)
      @:PROPAGATE_ERROR(errStatus)
      call getNodeName(value2, buffer)
      select case(char(buffer))
      case default
        call detailedError(child2, "Unknown method '" // char(buffer)&
            & // "' to generate electronegativities", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      case("paulingen")
        allocate(kENDefault(geo%nSpecies))
        kENDefault(:) = getElectronegativity(geo%speciesNames)
        call readSpeciesList(value2, geo%speciesNames, input%en, errStatus, kENDefault)
        @:PROPAGATE_ERROR(errStatus)
        deallocate(kENDefault)
      case("values")
        call readSpeciesList(value2, geo%speciesNames, input%en, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end select
      if (any(input%en <= 0.0_dp)) then
        call detailedError(value1, "Electronegativities are not defined for all species", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    else
      ! array is not used, but should still be populated with dummies
      input%en(:) = 0.0_dp
    end if

    allocate(input%covRad(geo%nSpecies))
    call getChildValue(value1, "Radii", value2, errStatus, "CovalentRadiiD3", child=child2)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value2, buffer)
    select case(char(buffer))
    case default
      call detailedError(child2, "Unknown method '"//char(buffer)//"' to generate radii", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case("covalentradiid3")
      allocate(kRadDefault(geo%nSpecies))
      kRadDefault(:) = getCovalentRadius(geo%speciesNames)
      call readSpeciesList(value2, geo%speciesNames, input%covRad, errStatus, kRadDefault)
      @:PROPAGATE_ERROR(errStatus)
      deallocate(kRadDefault)
    case("values")
      call readSpeciesList(value2, geo%speciesNames, input%covRad, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    if (any(input%covRad <= 0.0_dp)) then
      call detailedError(value1, "Covalent radii are not defined for all species", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine readCoordinationNumber


#:if WITH_MBD

  !> Reads in settings for Tkatchenko-Scheffler dispersion
  subroutine readDispTs(node, input, errStatus)

    !> data to parse
    type(fnode), pointer, intent(in) :: node

    !> control data coming back
    type(TDispMbdInp), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(string) :: buffer
    type(fnode), pointer :: child

    input%method = 'ts'
    call getChild(node, "EnergyAccuracy", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      call detailedWarning(child, "The energy accuracy setting will be ignored as it is not&
          & supported/need by libMBD any more")
    end if
    call getChild(node, "ForceAccuracy", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      call detailedWarning(child, "The force accuracy setting will be ignored as it is not&
          & supported/need by libMBD any more")
    end if
    call getChildValue(node, "Damping", input%ts_d, errStatus, default=(input%ts_d))
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "RangeSeparation", input%ts_sr, errStatus, default=(input%ts_sr))
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "ReferenceSet", buffer, errStatus, 'ts', child=child)
    @:PROPAGATE_ERROR(errStatus)
    input%vdw_params_kind = tolower(unquote(char(buffer)))
    call checkManyBodyDispRefName(input%vdw_params_kind, child, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "LogLevel", input%log_level, errStatus, default=(input%log_level))
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readDispTs


  !> Reads in many-body dispersion settings
  subroutine readDispMbd(node, input, errStatus)

    !> data to parse
    type(fnode), pointer, intent(in) :: node

    !> control data coming back
    type(TDispMbdInp), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(string) :: buffer
    type(fnode), pointer :: child

    input%method = 'mbd-rsscs'
    call getChildValue(node, "Beta", input%mbd_beta, errStatus, input%mbd_beta)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "NOmegaGrid", input%n_omega_grid, errStatus,&
        & default=(input%n_omega_grid))
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "KGrid", input%k_grid, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "KGridShift", input%k_grid_shift, errStatus,&
        & default=(input%k_grid_shift))
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "ReferenceSet", buffer, errStatus, 'ts', child=child)
    @:PROPAGATE_ERROR(errStatus)
    input%vdw_params_kind = tolower(unquote(char(buffer)))
    call checkManyBodyDispRefName(input%vdw_params_kind, child, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "LogLevel", input%log_level, errStatus, default=(input%log_level))
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readDispMbd


  !> Check the dispersion label matches allowed cases
  subroutine checkManyBodyDispRefName(name, node, errStatus)

    !> Label
    character(*), intent(in) :: name

    !> data tree for error usage
    type(fnode), pointer, intent(in) :: node

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    if (name /= 'ts' .and. name /= 'tssurf') then
      call detailedError(node, 'Invalid reference set name for TS/MBD-dispersion', errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine checkManyBodyDispRefName

#:endif

  !> reads in value of temperature for MD with sanity checking of the input
  subroutine readTemperature(node, ctrl, errStatus)

    !> data to parse
    type(fnode), pointer :: node

    !> control data coming back
    type(TControl), intent(inout) :: ctrl

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(string) :: modifier

    allocate(ctrl%tempSteps(1))
    allocate(ctrl%tempValues(1))
    allocate(ctrl%tempMethods(1))
    ctrl%tempMethods(1) = 1
    ctrl%tempSteps(1) = 1
    call getChildValue(node, "", ctrl%tempValues(1), errStatus, modifier=modifier)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), energyUnits, node, ctrl%tempValues(1), errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (ctrl%tempValues(1) < 0.0_dp) then
      call detailedError(node, "Negative temperature.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    if (ctrl%tempValues(1) < minTemp) then
      ctrl%tempValues(1) = minTemp
    end if

  end subroutine readTemperature


  !> reads a temperature profile for MD with sanity checking of the input
  subroutine readTemperatureProfile(node, modifier, ctrl, errStatus)

    !> parser node containing the relevant part of the user input
    type(fnode), pointer :: node

    !> unit modifier for the profile
    character(len=*), intent(in) :: modifier

    !> Control structure to populate
    type(TControl), intent(inout) :: ctrl

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(TListString) :: ls
    type(TListIntR1) :: li1
    type(TListRealR1) :: lr1
    character(len=20), allocatable :: tmpC1(:)
    integer :: ii
    logical :: success

    call init(ls)
    call init(li1)
    call init(lr1)
    call getChildValue(node, "", ls, 1, li1, 1, lr1, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (len(ls) < 1) then
      call detailedError(node, "At least one annealing step must be specified.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    allocate(tmpC1(len(ls)))
    allocate(ctrl%tempSteps(len(li1)))
    allocate(ctrl%tempValues(len(lr1)))
    call asArray(ls, tmpC1)
    call asVector(li1, ctrl%tempSteps)
    call asVector(lr1, ctrl%tempValues)
    call destruct(ls)
    call destruct(li1)
    call destruct(lr1)
    allocate(ctrl%tempMethods(size(tmpC1)))
    do ii = 1, size(tmpC1)
      call identifyTempProfile(ctrl%tempMethods(ii), tmpC1(ii), success)
      if (success) then
        cycle
      end if
      call detailedError(node, "Invalid annealing method name '" // trim(tmpC1(ii)) // "'.",&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end do

    if (any(ctrl%tempSteps < 0)) then
      call detailedError(node, "Step values must not be negative.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ii = sum(ctrl%tempSteps)
    if (ii < 1) then
      call detailedError(node, "Sum of steps in the profile must be greater than zero.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    ctrl%maxRun = ii - 1

    if (any(ctrl%tempValues < 0.0_dp)) then
      call detailedError(node, "Negative temperature.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call convertUnitHsd(modifier, energyUnits, node, ctrl%tempValues, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (any(ctrl%tempValues < minTemp)) then
      ctrl%tempValues = max(ctrl%tempValues, minTemp)
    end if
    deallocate(tmpC1)

  end subroutine readTemperatureProfile


  !> Reads the excited state data block
  subroutine readExcited(node, geo, ctrl, errStatus)

    !> Node to parse
    type(fnode), pointer :: node

    !> geometry object, which contains atomic species information
    type(TGeometry), intent(in) :: geo

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child
    type(fnode), pointer :: child2, child3
    type(fnode), pointer :: value
    type(string) :: buffer, modifier

    ! Linear response stuff
    call getChild(node, "Casida", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)

    if (associated(child)) then

      allocate(ctrl%lrespini)
      ctrl%lrespini%tPrintEigVecs = .false.

      if (ctrl%tSpin) then
        ctrl%lrespini%sym = ' '
      else
        call getChildValue(child, "Symmetry", buffer, errStatus, child=child2)
        @:PROPAGATE_ERROR(errStatus)
        select case (unquote(char(buffer)))
        case ("Singlet" , "singlet")
          ctrl%lrespini%sym = 'S'
        case ("Triplet" , "triplet")
          ctrl%lrespini%sym = 'T'
        case ("Both" , "both")
          ctrl%lrespini%sym = 'B'
        case default
          call detailedError(child2, "Invalid symmetry value '"  // char(buffer) // &
              & "' (must be 'Singlet', 'Triplet' or 'Both').", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end select
      end if

      call getChildValue(child, "NrOfExcitations", ctrl%lrespini%nexc, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChild(child, "StateOfInterest", child2, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (.not. associated(child2)) then
        ctrl%lrespini%nstat = 0
        call setChildValue(child, "StateOfInterest", 0, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      else
        call getChildValue(child2, "", buffer, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        if (tolower(unquote(char(buffer))) == "brightest") then
          if (ctrl%lrespini%sym /= "S" .or. ctrl%tSpin) then
            call detailedError(child2, "Brightest mode only allowed for spin unpolarised singlet&
                & excitations.", errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end if
          ctrl%lrespini%nstat = -1
        else
          call getChildValue(child2, "", ctrl%lrespini%nstat, errStatus)
          @:PROPAGATE_ERROR(errStatus)
          if (ctrl%lrespini%nstat > ctrl%lrespini%nexc) then
            call detailedError(child2, "Invalid value, must be within range of NrOfExcitations",&
                & errStatus)
            @:PROPAGATE_ERROR(errStatus)
          elseif (ctrl%lrespini%sym == "B" .and. ctrl%lrespini%nstat /= 0) then
            call detailedError(child2, "You cannot specify a particular excited state if symmetry&
                & is 'B'", errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end if
        end if
      end if

      call getChildValue(child, "EnergyWindow", ctrl%lrespini%energyWindow, errStatus, 0.0_dp,&
          & modifier=modifier, child=child2)
      @:PROPAGATE_ERROR(errStatus)
      ctrl%lrespini%tEnergyWindow = ctrl%lrespini%energyWindow /= 0.0_dp
      call convertUnitHsd(char(modifier), energyUnits, child2, ctrl%lrespini%energyWindow,&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child, "OscillatorWindow", ctrl%lrespini%oscillatorWindow, errStatus,&
          & 0.0_dp, modifier=modifier,  child=child2)
      @:PROPAGATE_ERROR(errStatus)
      ctrl%lrespini%tOscillatorWindow = ctrl%lrespini%oscillatorWindow /= 0.0_dp
      call convertUnitHsd(char(modifier), dipoleUnits, child2, ctrl%lrespini%oscillatorWindow,&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child, "CacheCharges", ctrl%lrespini%tCacheCharges, errStatus,&
          & default=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child, "WriteMulliken", ctrl%lrespini%tMulliken, errStatus,&
          & default=.false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child, "WriteCoefficients", ctrl%lrespini%tCoeffs, errStatus,&
          & default=.false.)
      @:PROPAGATE_ERROR(errStatus)
      ctrl%lrespini%tGrndState = .false.
      if (ctrl%lrespini%tCoeffs) then
        call getChildValue(child, "TotalStateCoeffs", ctrl%lrespini%tGrndState, errStatus, .false.)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getChildValue(child, "WriteEigenvectors", ctrl%lrespini%tPrintEigVecs, errStatus,&
          & .false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child, "WriteDensityMatrix", ctrl%lrespini%tWriteDensityMatrix, errStatus,&
          & .false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child, "WriteXplusY", ctrl%lrespini%tXplusY, errStatus, default=.false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child, "StateCouplings", ctrl%lrespini%indNACouplings, errStatus,&
          & default=[0, 0])
      @:PROPAGATE_ERROR(errStatus)
      if (all(ctrl%lrespini%indNACouplings == 0)) then
        ctrl%lrespini%tNaCoupling = .false.
      else
        ctrl%lrespini%tNaCoupling = .true.
      end if
      call getChildValue(child, "WriteSPTransitions", ctrl%lrespini%tSPTrans, errStatus,&
          & default=.false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child, "WriteTransitions", ctrl%lrespini%tTrans, errStatus,&
          & default=.false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child, "WriteTransitionDipole", ctrl%lrespini%tTradip, errStatus,&
          & default=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (allocated(ctrl%rangeSepInp)) then
        call getChildValue(child, "WriteTransitionCharges", ctrl%lrespini%tTransQ, errStatus,&
            & default=.false.)
        @:PROPAGATE_ERROR(errStatus)
      end if
      ctrl%lrespini%iLinRespSolver = linRespSolverTypes%None

      call renameChildren(child, "Diagonalizer", "Diagonaliser")
      call getChildValue(child, "Diagonaliser", child2, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getNodeName(child2, buffer)
      select case(char(buffer))
      case ("arpack")
        if (.not. withArpack) then
          call detailedError(child2, 'This DFTB+ binary has been compiled without support for&
              & linear response calculations using the ARPACK/ngARPACK library.', errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        call getChildValue(child2, "WriteStatusArnoldi", ctrl%lrespini%tArnoldi, errStatus,&
            & default=.false.)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child2, "TestArnoldi", ctrl%lrespini%tDiagnoseArnoldi, errStatus,&
            & default=.false.)
        @:PROPAGATE_ERROR(errStatus)
        ctrl%lrespini%iLinRespSolver = linRespSolverTypes%Arpack
      case ("stratmann")
        ctrl%lrespini%iLinRespSolver = linRespSolverTypes%Stratmann
        call getChildValue(child2, "SubSpaceFactor", ctrl%lrespini%subSpaceFactorStratmann,&
            & errStatus, 20)
        @:PROPAGATE_ERROR(errStatus)
      case default
        call detailedError(child2, "Invalid diagonaliser method '" // char(buffer) // "'",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end select

      call getChildValue(child, "OptimiserCI", child2, errStatus, "", child=child3,&
          & allowEmptyValue=.true.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child2)) then
        call getNodeName(child2, buffer)
        select case(char(buffer))
        case ("bearpark")
          ctrl%lrespini%isCIopt = .true.
          call getChildValue(child2, "EnergyShift", ctrl%lrespini%energyShiftCI,&
              & errStatus, modifier=modifier, default=0.0_dp)
          @:PROPAGATE_ERROR(errStatus)
          call convertUnitHsd(char(modifier), energyUnits, child, ctrl%lrespini%energyShiftCI,&
              & errStatus)
          @:PROPAGATE_ERROR(errStatus)
        case default
          call detailedError(child2, "Invalid optimiser method '" // char(buffer) // "'", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end select
      else
        ctrl%lrespini%isCIopt = .false.
      end if

      if (ctrl%tForces .or. ctrl%tPrintForces) then
        call getChildValue(child, "ExcitedStateForces", ctrl%tCasidaForces, errStatus,&
            & default=.true.)
        @:PROPAGATE_ERROR(errStatus)
      end if

    end if

    !pp-RPA
    call getChild(node, "PP-RPA", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)

    if (associated(child)) then

      allocate(ctrl%pprpa)

      if (ctrl%tSpin) then
        ctrl%pprpa%sym = ' '
      else
        call getChildValue(child, "Symmetry", buffer, errStatus, child=child2)
        @:PROPAGATE_ERROR(errStatus)
        select case (unquote(char(buffer)))
        case ("Singlet" , "singlet")
          ctrl%pprpa%sym = 'S'
        case ("Triplet" , "triplet")
          ctrl%pprpa%sym = 'T'
        case ("Both" , "both")
          ctrl%pprpa%sym = 'B'
        case default
          call detailedError(child2, "Invalid symmetry value '"  // char(buffer) // &
              & "' (must be 'Singlet', 'Triplet' or 'Both').", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end select
      end if

      call getChildValue(child, "NrOfExcitations", ctrl%pprpa%nexc, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(child, "HHubbard", value, errStatus, child=child2)
      @:PROPAGATE_ERROR(errStatus)
      allocate(ctrl%pprpa%hhubbard(geo%nSpecies))
      call readSpeciesList(child2, geo%speciesNames, ctrl%pprpa%hhubbard, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(child, "TammDancoff", ctrl%pprpa%tTDA, errStatus, default=.false.)
      @:PROPAGATE_ERROR(errStatus)

      call getChild(child, "NrOfVirtualStates", child2, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (.not. associated(child2)) then
        ctrl%pprpa%nvirtual = 0
        ctrl%pprpa%tConstVir = .false.
        call setChildValue(child, "NrOfVirtualStates", 0, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      else
        call getChildValue(child2, "", ctrl%pprpa%nvirtual, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        ctrl%pprpa%tConstVir = .true.
      end if

    end if

  end subroutine readExcited


  !> Reads the analysis block
#:if WITH_TRANSPORT
  subroutine readAnalysis(node, ctrl, geo, orb, transpar, tundos, errStatus)
#:else
  subroutine readAnalysis(node, ctrl, geo, errStatus)
#:endif

    !> Node to parse
    type(fnode), pointer :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

  #:if WITH_TRANSPORT
    !> Orbital
    type(TOrbitals), intent(in), allocatable :: orb

    !> Transport parameters
    type(TTransPar), intent(inout) :: transpar

    !> Tunneling and Dos parameters
    type(TNEGFTunDos), intent(inout) :: tundos
  #:endif

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: val, child, child2, child3
    type(fnodeList), pointer :: children
    integer, allocatable :: pTmpI1(:)
    type(string) :: buffer, modifier
    integer :: nReg, iReg
    character(lc) :: strTmp
    type(TListRealR1) :: lr1
    type(TListReal) :: lr
    logical :: tPipekDense
    logical :: tWriteBandDatDef, tHaveEigenDecomposition, tHaveDensityMatrix
    logical :: isEtaNeeded

    tHaveEigenDecomposition = .false.
    if (any(ctrl%solver%isolver == [electronicSolverTypes%qr,&
        & electronicSolverTypes%divideandconquer, electronicSolverTypes%relativelyrobust,&
        & electronicSolverTypes%elpa])) then
      tHaveEigenDecomposition = .true.
    end if
    tHaveDensityMatrix = ctrl%solver%isolver /= electronicSolverTypes%OnlyTransport

    if (tHaveEigenDecomposition) then

      call getChildValue(node, "ProjectStates", val, errStatus, "", child=child,&
          & allowEmptyValue=.true., list=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildren(child, "Region", children)
      nReg = getLength(children)
      ctrl%tProjEigenvecs = (nReg > 0)
      if (ctrl%tProjEigenvecs) then
        allocate(ctrl%tShellResInRegion(nReg))
        allocate(ctrl%tOrbResInRegion(nReg))
        allocate(ctrl%RegionLabel(nReg))
        call init(ctrl%iAtInRegion)
        do iReg = 1, nReg
          call getItem1(children, iReg, child2)
          call getChildValue(child2, "Atoms", buffer, errStatus, child=child3, multiple=.true.)
          @:PROPAGATE_ERROR(errStatus)
          call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species, pTmpI1,&
              & errStatus)
          @:PROPAGATE_ERROR(errStatus)
          call append(ctrl%iAtInRegion, pTmpI1)
          call getChildValue(child2, "ShellResolved", ctrl%tShellResInRegion(iReg), errStatus,&
              & .false., child=child3)
          @:PROPAGATE_ERROR(errStatus)
          if (ctrl%tShellResInRegion(iReg)) then
            if (.not. all(geo%species(pTmpI1) == geo%species(pTmpI1(1)))) then
              call detailedError(child3, "Shell resolved PDOS only allowed for &
                  &regions where all atoms belong to the same species", errStatus)
              @:PROPAGATE_ERROR(errStatus)
            end if
          end if
          call getChildValue(child2, "OrbitalResolved", ctrl%tOrbResInRegion(iReg), errStatus,&
              & .false., child=child3)
          @:PROPAGATE_ERROR(errStatus)
          if (ctrl%tOrbResInRegion(iReg)) then
            if (.not. all(geo%species(pTmpI1) == geo%species(pTmpI1(1)))) then
              call detailedError(child3, "Orbital resolved PDOS only allowed for &
                  &regions where all atoms belong to the same species", errStatus)
              @:PROPAGATE_ERROR(errStatus)
            end if
          end if
          deallocate(pTmpI1)
          write(strTmp, "('region',I0)") iReg
          call getChildValue(child2, "Label", buffer, errStatus, trim(strTmp))
          @:PROPAGATE_ERROR(errStatus)
          ctrl%RegionLabel(iReg) = unquote(char(buffer))
        end do
      end if
      call destroyNodeList(children)

      call renameChildren(node, "Localize", "Localise")
      call getChild(node, "Localise", val, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(val)) then
        ctrl%tLocalise = .true.
        call getChild(val, "PipekMezey", child2, errStatus, requested=.false.)
        @:PROPAGATE_ERROR(errStatus)
        if (associated(child2)) then
          allocate(ctrl%pipekMezeyInp)
          associate(inp => ctrl%pipekMezeyInp)
            call getChildValue(child2, "MaxIterations", inp%maxIter, errStatus, 100)
            @:PROPAGATE_ERROR(errStatus)
            tPipekDense = .true.
            if (.not. geo%tPeriodic) then
              call getChildValue(child2, "Dense", tPipekDense, errStatus, .false.)
              @:PROPAGATE_ERROR(errStatus)
              if (.not. tPipekDense) then
                call init(lr1)
                call getChild(child2, "SparseTolerances", child3, errStatus, requested=.false.)
                @:PROPAGATE_ERROR(errStatus)
                if (associated(child3)) then
                  call getChildValue(child3, "", 1, lr1, errStatus)
                  @:PROPAGATE_ERROR(errStatus)
                  if (len(lr1) < 1) then
                    call detailedError(child2, "Missing values of tolerances.", errStatus)
                    @:PROPAGATE_ERROR(errStatus)
                  end if
                  allocate(inp%sparseTols(len(lr1)))
                  call asVector(lr1, inp%sparseTols)
                else
                  allocate(inp%sparseTols(4))
                  inp%sparseTols = [0.1_dp, 0.01_dp, 1.0E-6_dp, 1.0E-12_dp]
                  call setChildValue(child2, "SparseTolerances", inp%sparseTols, errStatus)
                  @:PROPAGATE_ERROR(errStatus)
                end if
                call destruct(lr1)
              end if
            end if
            if (tPipekDense) then
              call getChildValue(child2, "Tolerance", inp%tolerance, errStatus, 1.0E-4_dp)
              @:PROPAGATE_ERROR(errStatus)
            end if
          end associate
        else
          call detailedError(val, "No localisation method chosen", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      end if

      call getChildValue(node, "WriteEigenvectors", ctrl%tPrintEigVecs, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)

    #:if WITH_SOCKETS
      tWriteBandDatDef = .not. allocated(ctrl%socketInput)
    #:else
      tWriteBandDatDef = .true.
    #:endif

      call getChildValue(node, "WriteBandOut", ctrl%tWriteBandDat, errStatus, tWriteBandDatDef)
      @:PROPAGATE_ERROR(errStatus)

      call getChild(node, "Polarisability", child, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChild(node, "ResponseKernel", child2, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child) .or. associated(child2)) then
        allocate(ctrl%perturbInp)
      end if

      ! electric field polarisability of system
      call getChild(node, "Polarisability", child, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child)) then
        ctrl%perturbInp%isEPerturb = .true.
        call freqRanges(child, ctrl%perturbInp%dynEFreq, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

      call getChild(node, "ResponseKernel", child, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(child)) then
        ctrl%perturbInp%isRespKernelPert = .true.
        if (ctrl%tSCC) then
          call getChildValue(child, "RPA", ctrl%perturbInp%isRespKernelRPA, errStatus, .false.)
          @:PROPAGATE_ERROR(errStatus)
        else
          ctrl%perturbInp%isRespKernelRPA = .true.
        end if
        call freqRanges(child, ctrl%perturbInp%dynKernelFreq, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

      if (allocated(ctrl%perturbInp)) then
        call getChildValue(node, "PertubDegenTol", ctrl%perturbInp%tolDegenDFTBPT, errStatus,&
            & 128.0_dp, child=child)
        @:PROPAGATE_ERROR(errStatus)
        if (ctrl%perturbInp%tolDegenDFTBPT < 1.0_dp) then
          call detailedError(child, "Perturbation degeneracy tolerance must be above 1x", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        ctrl%perturbInp%tolDegenDFTBPT = ctrl%perturbInp%tolDegenDFTBPT * epsilon(0.0_dp)
        isEtaNeeded = .false.
        if (allocated(ctrl%perturbInp%dynEFreq)) then
          if (any(ctrl%perturbInp%dynEFreq /= 0.0_dp)) then
            isEtaNeeded = .true.
          end if
        end if
        if (allocated(ctrl%perturbInp%dynKernelFreq)) then
          if (any(ctrl%perturbInp%dynKernelFreq /= 0.0_dp)) then
            isEtaNeeded = .true.
          end if
        end if
        if (isEtaNeeded) then
          allocate(ctrl%perturbInp%etaFreq)
          call getChildValue(node, "PerturbEta", ctrl%perturbInp%etaFreq, errStatus, 1.0E-8_dp,&
              & child=child)
          @:PROPAGATE_ERROR(errStatus)
          if (ctrl%perturbInp%etaFreq < epsilon(0.0_dp)) then
            call detailedError(child, "Imaginary constant for finite frequency perturbation too&
                & small", errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end if
        end if
      end if

      if (allocated(ctrl%perturbInp)) then
        call maxSelfConsIterations(node, ctrl, "MaxPerturbIter", ctrl%perturbInp%maxPerturbIter,&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
        if (ctrl%tScc) then
          call getChildValue(node, "PerturbSccTol", ctrl%perturbInp%perturbSccTol, errStatus,&
              & 1.0e-5_dp)
          @:PROPAGATE_ERROR(errStatus)
          ! self consistency required, or not, to proceed with perturbation
          call getChildValue(node, "ConvergedPerturb", ctrl%perturbInp%isPerturbConvRequired,&
              & errStatus, .true.)
          @:PROPAGATE_ERROR(errStatus)
        end if
      end if

    end if

    if (tHaveDensityMatrix) then

      ! Is this compatible with Poisson solver use?
      call readElectrostaticPotential(node, geo, ctrl, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(node, "MullikenAnalysis", ctrl%tPrintMulliken, errStatus, .true.)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%tPrintMulliken) then
        call getChildValue(node, "WriteNetCharges", ctrl%tPrintNetAtomCharges, errStatus,&
            & default=.false.)
        @:PROPAGATE_ERROR(errStatus)
        if (ctrl%tPrintNetAtomCharges) then
          ctrl%tNetAtomCharges = .true.
        end if
        call getChild(node, "CM5", child, errStatus, requested=.false.)
        @:PROPAGATE_ERROR(errStatus)
        if (associated(child)) then
          allocate(ctrl%cm5Input)
          call readCM5(child, ctrl%cm5Input, geo, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      end if
      call getChildValue(node, "AtomResolvedEnergies", ctrl%tAtomicEnergy, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)

      if (allocated(ctrl%solvInp)) then
        call getChildValue(node, "writeCosmoFile", ctrl%tWriteCosmoFile, errStatus,&
            & allocated(ctrl%solvInp%cosmoInp), child=child)
        @:PROPAGATE_ERROR(errStatus)
        if (ctrl%tWriteCosmoFile .and. .not.allocated(ctrl%solvInp%cosmoInp)) then
          call detailedError(child, "Cosmo file can only be written for Cosmo calculations",&
              & errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      end if

      call getChildValue(node, "PrintForces", ctrl%tPrintForces, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)

    else

      ctrl%tPrintMulliken = .false.
      ctrl%tAtomicEnergy = .false.
      ctrl%tPrintForces = .false.

    end if


  #:if WITH_TRANSPORT
    call getChild(node, "TunnelingAndDOS", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      if (.not.transpar%defined) then
        @:RAISE_ERROR(errStatus, -1, "Block TunnelingAndDos requires Transport block.")
      end if
      if (.not.transpar%taskUpload) then
        @:RAISE_ERROR(errStatus, -1, "Block TunnelingAndDos not compatible with&
            & task=contactHamiltonian")
      end if
      if (.not. allocated(orb)) then
        @:RAISE_ERROR(errStatus, -1, "Orbital information from SK-files missing (xTB Hamiltonian&
            & not compatible with transport yet)")
      end if
      call readTunAndDos(child, orb, geo, tundos, transpar, ctrl%tempElec, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else
      if (ctrl%solver%isolver == electronicSolverTypes%OnlyTransport) then
        call detailedError(node, "The TransportOnly solver requires a TunnelingAndDos block to be&
            & present.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    endif
  #:endif

  end subroutine readAnalysis


  !> Frequency ranges for response calculations
  subroutine freqRanges(node, frequencies, errStatus)

    !> Node to parse
    type(fnode), pointer :: node

    !> Frequencies, 0 being static
    real(dp), allocatable, intent(inout) :: frequencies(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(TListReal) :: lr
    type(fnode), pointer :: child, child2
    type(string) :: modifier
    integer :: nFreq, iFreq, jFreq
    real(dp) :: tmp3R(3)
    logical :: isStatic

    call getChildValue(node, "Static", isStatic, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)
    if (isStatic) then
      call growFreqArray(frequencies, 1)
      ! should already be zero, but just in case:
      frequencies(:) = 0.0_dp
    end if

    call getChild(node, "Frequencies", child, errStatus, modifier=modifier, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      call init(lr)
      call getChildValue(child, "", lr, errStatus, child=child2, modifier=modifier)
      @:PROPAGATE_ERROR(errStatus)
      nFreq = len(lr)
      if (nFreq > 0) then
        if (allocated(frequencies)) then
          iFreq = size(frequencies)
        else
          iFreq = 0
        end if
        call growFreqArray(frequencies, nFreq)
        call asArray(lr, frequencies(iFreq+1:iFreq+nFreq))
        call convertUnitHsd(char(modifier),freqUnits, child, frequencies(iFreq+1:iFreq+nFreq),&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call destruct(lr)
      if (any(frequencies < 0.0_dp)) then
        call detailedError(child2, "Negative driving frequency requested", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end if

    call getChild(node, "FrequencyRange", child, errStatus, modifier=modifier, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      call init(lr)
      call getChildValue(child, "", lr, errStatus, child=child2, modifier=modifier)
      @:PROPAGATE_ERROR(errStatus)
      if (len(lr) == 3) then
        call asArray(lr, tmp3R)
        call convertUnitHsd(char(modifier), freqUnits, child, tmp3R, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        if (any(tmp3R(:2) < 0.0_dp)) then
          call detailedError(child, "Negative values in dynamic frequency range.", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        if (abs(tmp3R(3)) <= epsilon(0.0_dp)) then
          call detailedError(child, "Increase step size in dynamic frequency range.", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        ! how many frequencies in the specified range?
        nFreq = max(int((tmp3R(2)-tmp3R(1))/tmp3R(3))+1,0)
        if (allocated(frequencies)) then
          iFreq = size(frequencies)
        else
          iFreq = 0
        end if
        call growFreqArray(frequencies, nFreq)
        do jFreq = 1, nFreq
          frequencies(iFreq+jFreq) = tmp3R(1) + (jFreq-1) * tmp3R(3)
        end do
      else
        call detailedError(child,"Malformed frequency range.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call destruct(lr)
      if (any(frequencies < 0.0_dp)) then
        call detailedError(child2, "Negative driving frequency requested", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end if

  end subroutine freqRanges


  !> Resize array, retaining values at start
  subroutine growFreqArray(freq, nFreq)

    !> Array to expand
    real(dp), allocatable, intent(inout) :: freq(:)

    !> Number of extra elements
    integer, intent(in) :: nFreq

    real(dp), allocatable :: tmpFreq(:)
    integer :: nElem

    if (allocated(freq)) then
      nElem = size(freq)
      call move_alloc(freq, tmpFreq)
    else
      nElem =0
    end if
    allocate(freq(nElem + nFreq))
    if (nElem > 0) then
      freq(:nElem) = tmpFreq
    end if
    freq(nElem+1:) = 0.0_dp

  end subroutine growFreqArray


  !> Read in settings that are influenced by those read from Options{} but belong in Analysis{}
  subroutine readLaterAnalysis(node, ctrl, errStatus)

    !> Node to parse
    type(fnode), pointer :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> Error status
    type(TStatus), intent(inout) :: errStatus


    logical :: tPrintEigVecs

    tPrintEigVecs = ctrl%tPrintEigVecs
    if (allocated(ctrl%lrespini)) tPrintEigvecs = tPrintEigvecs .or. ctrl%lrespini%tPrintEigVecs
    if (tPrintEigVecs) then
      call getChildValue(node, "EigenvectorsAsText", ctrl%tPrintEigVecsTxt, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine readLaterAnalysis


  !> Read in hamiltonian settings that are influenced by those read from REKS{}, electronDynamics{}
  subroutine readLaterHamiltonian(hamNode, ctrl, driverNode, geo, errStatus)

    !> Hamiltonian node to parse
    type(fnode), pointer :: hamNode

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> Geometry driver node to parse
    type(fnode), pointer :: driverNode

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, value2, child, child2
    type(string) :: buffer, buffer2
    type(TListRealR1) :: lr1

    if (ctrl%reksInp%reksAlg == reksTypes%noReks) then

      if (ctrl%tSCC) then

        call getChildValue(hamNode, "Mixer", value1, errStatus, "Broyden", child=child)
        @:PROPAGATE_ERROR(errStatus)
        call getNodeName(value1, buffer)
        select case(char(buffer))

        case ("broyden")

          ctrl%iMixSwitch = mixerTypes%broyden
          call getChildValue(value1, "MixingParameter", ctrl%almix, errStatus, 0.2_dp)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(value1, "InverseJacobiWeight", ctrl%broydenOmega0, errStatus, 0.01_dp)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(value1, "MinimalWeight", ctrl%broydenMinWeight, errStatus, 1.0_dp)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(value1, "MaximalWeight", ctrl%broydenMaxWeight, errStatus, 1.0e5_dp)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(value1, "WeightFactor", ctrl%broydenWeightFac, errStatus, 1.0e-2_dp)
          @:PROPAGATE_ERROR(errStatus)

        case ("anderson")

          ctrl%iMixSwitch = mixerTypes%anderson
          call getChildValue(value1, "MixingParameter", ctrl%almix, errStatus, 0.05_dp)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(value1, "Generations", ctrl%iGenerations, errStatus, 4)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(value1, "InitMixingParameter", ctrl%andersonInitMixing, errStatus,&
              & 0.01_dp)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(value1, "DynMixingParameters", value2, errStatus, "", child=child,&
              & allowEmptyValue=.true.)
          @:PROPAGATE_ERROR(errStatus)
          call getNodeName2(value2, buffer2)
          if (char(buffer2) == "") then
            ctrl%andersonNrDynMix = 0
          else
            call init(lr1)
            call getChildValue(child, "", 2, lr1, errStatus, child=child2)
            @:PROPAGATE_ERROR(errStatus)
            if (len(lr1) < 1) then
              call detailedError(child2, "At least one dynamic mixing parameter must be defined.",&
                  & errStatus)
              @:PROPAGATE_ERROR(errStatus)
            end if
            ctrl%andersonNrDynMix = len(lr1)
            allocate(ctrl%andersonDynMixParams(2, ctrl%andersonNrDynMix))
            call asArray(lr1, ctrl%andersonDynMixParams)
            call destruct(lr1)
          end if
          call getChildValue(value1, "DiagonalRescaling", ctrl%andersonOmega0, errStatus, 1.0e-2_dp)
          @:PROPAGATE_ERROR(errStatus)

        case ("simple")

          ctrl%iMixSwitch = mixerTypes%simple
          call getChildValue(value1, "MixingParameter", ctrl%almix, errStatus, 0.05_dp)
          @:PROPAGATE_ERROR(errStatus)

        case("diis")

          ctrl%iMixSwitch = mixerTypes%diis
          call getChildValue(value1, "InitMixingParameter", ctrl%almix, errStatus, 0.2_dp)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(value1, "Generations", ctrl%iGenerations, errStatus, 6)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(value1, "UseFromStart", ctrl%tFromStart, errStatus, .true.)
          @:PROPAGATE_ERROR(errStatus)

        case default

          call getNodeHSDName(value1, buffer)
          call detailedError(child, "Invalid mixer '" // char(buffer) // "'", errStatus)
          @:PROPAGATE_ERROR(errStatus)

        end select

      end if

      if (ctrl%tMD) then
        if (ctrl%iThermostat /= 0) then
          call getChildValue(driverNode, "Thermostat", child, errStatus, child=child2)
          @:PROPAGATE_ERROR(errStatus)
          if (ctrl%reksInp%reksAlg == reksTypes%noReks) then
            call getChildValue(child, "AdaptFillingTemp", ctrl%tSetFillingTemp, errStatus, .false.)
            @:PROPAGATE_ERROR(errStatus)
          end if
        end if
      end if

    end if

    hamNeedsT: if (ctrl%reksInp%reksAlg == reksTypes%noReks) then

      if (allocated(ctrl%elecDynInp)) then
        if (ctrl%elecDynInp%tReadRestart .and. .not.ctrl%elecDynInp%tPopulations) then
          exit hamNeedsT
        end if
      end if

      if (ctrl%solver%isolver /= electronicSolverTypes%GF) then
        call readElectronicFilling(hamNode, ctrl, geo, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

    end if hamNeedsT

  end subroutine readLaterHamiltonian


  !> Parses for electronic filling temperature (should only read if not either REKS or electron
  !> dynamics from a supplied density matrix)
  subroutine readElectronicFilling(hamNode, ctrl, geo, errStatus)

    !> Relevant node in input tree
    type(fnode), pointer :: hamNode

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to test for periodicity
    type(TGeometry), intent(in) :: geo

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    select case(ctrl%hamiltonian)
    case(hamiltonianTypes%xtb)
      call readFilling(hamNode, ctrl, geo, 300.0_dp*Boltzmann, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case(hamiltonianTypes%dftb)
      call readFilling(hamNode, ctrl, geo, 0.0_dp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readElectronicFilling


  !> Reads W values if required by settings in the Hamiltonian or the excited state
  subroutine readSpinConstants(hamNode, geo, orb, ctrl, errStatus)

    !> node for Hamiltonian data
    type(fnode), pointer :: hamNode

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> control structure
    type(TControl), intent(inout) :: ctrl

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child
    logical :: tLRNeedsSpinConstants, tShellResolvedW
    integer :: iSp1

    tLRNeedsSpinConstants = .false.

    if (allocated(ctrl%lrespini)) then
      select case (ctrl%lrespini%sym)
      case ("T", "B", " ")
        tLRNeedsSpinConstants = .true.
      case ("S")
        tLRNeedsSpinConstants = .false.
      case default
      end select
    end if

    if (tLRNeedsSpinConstants .or. ctrl%tSpin .or. &
        & ctrl%reksInp%reksAlg /= reksTypes%noReks) then
      allocate(ctrl%spinW(orb%mShell, orb%mShell, geo%nSpecies))
      ctrl%spinW(:,:,:) = 0.0_dp

      call getChild(hamNode, "SpinConstants", child, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%hamiltonian == hamiltonianTypes%xtb) then
        call getChildValue(child, "ShellResolvedSpin", tShellResolvedW, errStatus, .true.)
        @:PROPAGATE_ERROR(errStatus)
      else
        if (.not.ctrl%tShellResolved) then
          call getChildValue(child, "ShellResolvedSpin", tShellResolvedW, errStatus, .false.)
          @:PROPAGATE_ERROR(errStatus)
        else
          tShellResolvedW = .true.
        end if
      end if

      if (tShellResolvedW) then
        ! potentially unique values for each shell
        do iSp1 = 1, geo%nSpecies
          call getChildValue(child, geo%speciesNames(iSp1),&
              & ctrl%spinW(:orb%nShell(iSp1), :orb%nShell(iSp1), iSp1), errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end do
      else
        ! only one value per atom
        do iSp1 = 1, geo%nSpecies
          call getChildValue(child, geo%speciesNames(iSp1), ctrl%spinW(1, 1, iSp1), errStatus)
          @:PROPAGATE_ERROR(errStatus)
          ctrl%spinW(:orb%nShell(iSp1), :orb%nShell(iSp1), iSp1) =&
              & ctrl%spinW(1, 1, iSp1)
        end do
      end if
    end if

  end subroutine readSpinConstants


  !> Reads customised Hubbard U values that over-ride the SK file values
  subroutine readCustomisedHubbards(node, geo, orb, tShellResolvedScc, hubbU, errStatus)

    !> input data to parse
    type(fnode), pointer, intent(in) :: node

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> is this a shell resolved calculation, or only one U value per atom
    logical, intent(in) :: tShellResolvedScc

    !> hubbard U values on exit
    real(dp), allocatable, intent(out) :: hubbU(:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child, child2
    integer :: iSp1

    call renameChildren(node, "CustomizedHubbards", "CustomisedHubbards")
    call getChild(node, "CustomisedHubbards", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      allocate(hubbU(orb%mShell, geo%nSpecies))
      hubbU(:,:) = 0.0_dp
      do iSp1 = 1, geo%nSpecies
        call getChild(child, geo%speciesNames(iSp1), child2, errStatus, requested=.false.)
        @:PROPAGATE_ERROR(errStatus)
        if (.not. associated(child2)) then
          cycle
        end if
        if (tShellResolvedScc) then
          call getChildValue(child2, "", hubbU(:orb%nShell(iSp1), iSp1), errStatus)
          @:PROPAGATE_ERROR(errStatus)
        else
          call getChildValue(child2, "", hubbU(1, iSp1), errStatus)
          @:PROPAGATE_ERROR(errStatus)
          hubbU(:orb%nShell(iSp1), iSp1) = hubbU(1, iSp1)
        end if
      end do
    end if

  end subroutine readCustomisedHubbards


  !> Reads the electron dynamics block
  subroutine readElecDynamics(node, input, geom, masses, errStatus)

    !> input data to parse
    type(fnode), pointer :: node

    !> ElecDynamicsInp instance
    type(TElecDynamicsInp), intent(inout) :: input

    !> geometry of the system
    type(TGeometry), intent(in) :: geom

    !> masses to be returned
    real(dp), allocatable, intent(inout) :: masses(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child
    type(string) :: buffer, buffer2, modifier
    logical :: ppRangeInvalid, tNeedFieldStrength
    real (dp) :: defPpRange(2)
    logical :: defaultWrite

  #:if WITH_MPI
    if (associated(node)) then
      call detailedError(node, 'This DFTB+ binary has been compiled with MPI settings and &
          & electron dynamics are not currently available for distributed parallel calculations.',&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
  #:endif

    call getChildValue(node, "Steps", input%steps, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "TimeStep", input%dt, errStatus, modifier=modifier, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), timeUnits, child, input%dt, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(node, "Populations", input%tPopulations, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "WriteFrequency", input%writeFreq, errStatus, 50)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "Restart", input%tReadRestart, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    if (input%tReadRestart) then
      call getChildValue(node, "RestartFromAscii", input%tReadRestartAscii, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
    end if
    call getChildValue(node, "WriteRestart", input%tWriteRestart, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)
    if (input%tWriteRestart) then
      call getChildValue(node, "WriteAsciiRestart", input%tWriteRestartAscii, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
    end if
    call getChildValue(node, "RestartFrequency", input%restartFreq, errStatus,&
        & max(input%Steps / 10, 1))
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "Forces", input%tForces, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "WriteBondEnergy", input%tBondE, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "WriteBondPopulation", input%tBondP, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "WriteAtomicEnergies", input%tWriteAtomEnergies, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "Pump", input%tPump, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "FillingsFromFile", input%tFillingsFromFile, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)

    if (input%tPump) then
      call getChildValue(node, "PumpProbeFrames", input%tdPPFrames, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      defPpRange = [0.0_dp, input%steps * input%dt]
      call getChildValue(node, "PumpProbeRange", input%tdPpRange, errStatus, defPprange,&
          & modifier=modifier, child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), timeUnits, child, input%tdPpRange, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      ppRangeInvalid = (input%tdPpRange(2) <= input%tdPpRange(1))&
          & .or. (input%tdPprange(1) < defPpRange(1))&
          & .or. (input%tdPpRange(2) > defPpRange(2))
      if (ppRangeInvalid) then
        call detailederror(child, "Wrong definition of PumpProbeRange, either incorrect order&
            & or outside of simulation time range", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end if

    call getChildValue(node, "Probe", input%tProbe, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    if (input%tPump .and. input%tProbe) then
      call detailedError(child, "Pump and probe cannot be simultaneously true.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getChildValue(node, "EulerFrequency", input%eulerFreq, errStatus, 0)
    @:PROPAGATE_ERROR(errStatus)
    if ((input%eulerFreq < 50) .and. (input%eulerFreq > 0)) then
      call detailedError(child, "Wrong number of Euler steps, should be above 50", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    if (input%eulerFreq >= 50) then
      input%tEulers = .true.
    else
      input%tEulers = .false.
    end if

    ! assume this is required (needed for most perturbations, but not none)
    tNeedFieldStrength = .true.

    defaultWrite = .true.

    !! Different perturbation types
    call getChildValue(node, "Perturbation", value1, errStatus, "None", child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    select case(char(buffer))

    case ("kick")
      input%pertType = pertTypes%kick
      call renameChildren(value1, "PolarizationDirection", "PolarisationDirection")
      call getChildValue(value1, "PolarisationDirection", buffer2, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call directionConversion(unquote(char(buffer2)), value1, input%polDir, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(value1, "SpinType", buffer2, errStatus, "Singlet")
      @:PROPAGATE_ERROR(errStatus)
      select case(unquote(char(buffer2)))
      case ("singlet", "Singlet")
        input%spType = tdSpinTypes%singlet
      case ("triplet", "Triplet")
        input%spType = tdSpinTypes%triplet
      case default
        call detailedError(value1, "Unknown spectrum spin type " // char(buffer2), errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end select

      defaultWrite = .false.

    case ("laser")
      input%pertType = pertTypes%laser
      call renameChildren(value1, "PolarizationDirection", "PolarisationDirection")
      call getChildValue(value1, "PolarisationDirection", input%reFieldPolVec, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call renameChildren(value1, "ImagPolarizationDirection", "ImagPolarisationDirection")
      call getChildValue(value1, "ImagPolarisationDirection", input%imFieldPolVec, errStatus,&
          & [0.0_dp, 0.0_dp, 0.0_dp])
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "LaserEnergy", input%omega, errStatus, modifier=modifier,&
          & child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, child, input%omega, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "Phase", input%phase, errStatus, 0.0_dp, modifier=modifier,&
          & child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), angularUnits, child, input%phase, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "ExcitedAtoms", buffer, errStatus, "1:-1", child=child,&
          & multiple=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getSelectedAtomIndices(child, char(buffer), geom%speciesNames, geom%species,&
          & input%indExcitedAtom, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      input%nExcitedAtom = size(input%indExcitedAtom)
      if (input%nExcitedAtom == 0) then
        @:RAISE_ERROR(errStatus, -1, "No atoms specified for laser excitation.")
      end if

      defaultWrite = .true.

    case ("kickandlaser")
      input%pertType = pertTypes%kickAndLaser
      call getChildValue(value1, "KickPolDir", buffer2, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call directionConversion(unquote(char(buffer2)), value1, input%polDir, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "SpinType", input%spType, errStatus, tdSpinTypes%singlet)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "LaserPolDir", input%reFieldPolVec, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "LaserImagPolDir", input%imFieldPolVec, errStatus,&
          & [0.0_dp, 0.0_dp, 0.0_dp])
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "LaserEnergy", input%omega, errStatus, modifier=modifier,&
          & child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, child, input%omega, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "Phase", input%phase, errStatus, 0.0_dp, modifier=modifier,&
          & child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), angularUnits, child, input%phase, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(value1, "LaserStrength", input%tdLaserField, errStatus, modifier=modifier,&
          & child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), EFieldUnits, child, input%tdLaserField, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(value1, "ExcitedAtoms", buffer, errStatus, "1:-1", child=child,&
          & multiple=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getSelectedAtomIndices(child, char(buffer), geom%speciesNames, geom%species,&
          & input%indExcitedAtom, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      input%nExcitedAtom = size(input%indExcitedAtom)
      if (input%nExcitedAtom == 0) then
        @:RAISE_ERROR(errStatus, -1, "No atoms specified for laser excitation.")
      end if

      defaultWrite = .false.

    case ("none")
      input%pertType = pertTypes%noTDPert
      tNeedFieldStrength = .false.

      defaultWrite = .true.

    case default
      call detailedError(child, "Unknown perturbation type " // char(buffer), errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    if (tNeedFieldStrength) then
      call getChildValue(node, "FieldStrength", input%tdfield, errStatus, modifier=modifier,&
          & child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), EFieldUnits, child, input%tdfield, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getChildValue(node, "WriteEnergyAndCharges", input%tdWriteExtras, errStatus, defaultWrite)
    @:PROPAGATE_ERROR(errStatus)

    !! Different envelope functions
    call getChildValue(node, "EnvelopeShape", value1, errStatus, "Constant")
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value1, buffer)
    select case(char(buffer))

    case("constant")
      input%envType = envTypes%constant

    case("gaussian")
      input%envType = envTypes%gaussian
      call getChildValue(value1, "Time0", input%time0, errStatus, 0.0_dp, modifier=modifier,&
          & child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), timeUnits, child, input%Time0, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(value1, "Time1", input%time1, errStatus, modifier=modifier, child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), timeUnits, child, input%Time1, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case("sin2")
      input%envType = envTypes%sin2
      call getChildValue(value1, "Time0", input%time0, errStatus, 0.0_dp, modifier=modifier,&
          & child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), timeUnits, child, input%Time0, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(value1, "Time1", input%time1, errStatus, modifier=modifier, child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), timeUnits, child, input%Time1, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case("fromfile")
      input%envType = envTypes%fromFile
      call getChildValue(value1, "Time0", input%time0, errStatus, 0.0_dp, modifier=modifier,&
          & child=child)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), timeUnits, child, input%Time0, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case default
      call detailedError(value1, "Unknown envelope shape " // char(buffer), errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    !! Non-adiabatic molecular dynamics
    call getChildValue(node, "IonDynamics", input%tIons, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    if (input%tIons) then
      call getChildValue(node, "MovedAtoms", buffer, errStatus, "1:-1", child=child,&
          & multiple=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getSelectedAtomIndices(child, char(buffer), geom%speciesNames, geom%species,&
          & input%indMovedAtom, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      input%nMovedAtom = size(input%indMovedAtom)
      call readInitialVelocitiesNAMD(node, input, geom%nAtom, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (input%tReadMDVelocities) then
        ! without a thermostat, if we know the initial velocities, we do not need a temperature, so
        ! just set it to something 'safe'
        input%tempAtom = minTemp
      else
        if (.not. input%tReadRestart) then
          ! previously lower limit was minTemp:
          call readMDInitTemp(node, input%tempAtom, 0.0_dp, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        call getInputMasses(node, geom, masses, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end if

  end subroutine readElecDynamics


  !> Read in initial ion temperature for simple MD
  subroutine readMDInitTemp(node, tempAtom, minimumTemp, errStatus)

    !> input data to parse
    type(fnode), pointer :: node

    !> Ionic temperature
    real(dp), intent(out) :: tempAtom

    !> Lowest possible ion temperature
    real(dp), intent(in) :: minimumTemp

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child
    type(string) :: modifier

    call getChildValue(node, "InitialTemperature", tempAtom, errStatus, modifier=modifier,&
        & child=child)
    @:PROPAGATE_ERROR(errStatus)
    if (tempAtom < 0.0_dp) then
      call detailedError(node, "Negative temperature", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    call convertUnitHsd(char(modifier), energyUnits, node, tempAtom, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    tempAtom = max(tempAtom, minimumTemp)

  end subroutine readMDInitTemp


  !> Converts direction label text string into corresponding numerical value
  subroutine directionConversion(direction, node, iX, errStatus)

    !> Direction label
    character(*), intent(in) :: direction

    !> Input tree for error return
    type(fnode), pointer :: node

    !> Direction indicator (1 - 4) for (x,y,z,all)
    integer, intent(out) :: iX

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    select case(trim(direction))
    case ("x", "X")
      iX = 1
    case ("y", "Y")
      iX = 2
    case ("z", "Z")
      iX = 3
    case ("all", "All", "ALL")
      iX = 4
    case default
      call detailedError(node, "Wrongly specified polarisation direction " // trim(direction)&
          & // ". Must be x, y, z or all.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine directionConversion


  !> Reads MD velocities
  subroutine readInitialVelocitiesNAMD(node, input, nAtom, errStatus)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> ElecDynamicsInp object structure to be filled
    type(TElecDynamicsInp), intent(inout) :: input

    !> Total number of all atoms
    integer, intent(in) :: nAtom

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child
    type(string) :: buffer, modifier
    type(TListRealR1) :: realBuffer
    integer :: nVelocities
    real(dp), pointer :: tmpVelocities(:,:)

    call getChildValue(node, "Velocities", value1, errStatus, "", child=child,&
        & modifier=modifier, allowEmptyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName2(value1, buffer)
    if (char(buffer) == "") then
      input%tReadMDVelocities = .false.
    else
      call init(realBuffer)
      call getChildValue(child, "", 3, realBuffer, errStatus, modifier=modifier)
      @:PROPAGATE_ERROR(errStatus)
      nVelocities = len(realBuffer)
      if (nVelocities /= nAtom) then
        call detailedError(node, "Incorrect number of specified velocities: "&
            & // i2c(3*nVelocities) // " supplied, "&
            & // i2c(3*nAtom) // " required.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      allocate(tmpVelocities(3, nVelocities))
      call asArray(realBuffer, tmpVelocities)
      if (len(modifier) > 0) then
        call convertUnitHsd(char(modifier), VelocityUnits, child, tmpVelocities, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call destruct(realBuffer)
      allocate(input%initialVelocities(3, input%nMovedAtom))
      input%initialVelocities(:,:) = tmpVelocities(:, input%indMovedAtom(:))
      input%tReadMDVelocities = .true.
    end if

  end subroutine readInitialVelocitiesNAMD


#:if WITH_TRANSPORT
  !> Read geometry information for transport calculation
  subroutine readTransportGeometry(root, geom, transpar, errStatus)

    !> Root node containing the current block
    type(fnode), pointer :: root

    !> geometry of the system, which may be modified for some types of calculation
    type(TGeometry), intent(inout) :: geom

    !> Parameters of the transport calculation
    type(TTransPar), intent(inout) :: transpar

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: pDevice, pTask, pTaskType
    type(string) :: buffer, modifier
    type(fnode), pointer :: pTmp, field
    type(fnodeList), pointer :: pNodeList
    integer :: contact
    real(dp) :: lateralContactSeparation
    logical, allocatable :: atomInRegion(:)
    integer :: ii
    character(lc) :: strTmp

    transpar%defined = .true.
    transpar%tPeriodic1D = .not. geom%tPeriodic

    !! Note: we parse first the task because we need to know it to define the
    !! mandatory contact entries. On the other hand we need to wait that
    !! contacts are parsed to resolve the name of the contact for task =
    !! contacthamiltonian
    call getChildValue(root, "Task", pTaskType, errStatus, child=pTask, default='uploadcontacts')
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(pTaskType, buffer)

    call getChild(root, "Device", pDevice, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pDevice, "AtomRange", transpar%idxdevice, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChild(pDevice, "FirstLayerAtoms", pTmp, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    call readFirstLayerAtoms(pTmp, transpar%PL, transpar%nPLs, transpar%idxdevice, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (.not.associated(pTmp)) then
      call setChildValue(pDevice, "FirstLayerAtoms", transpar%PL, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getChildren(root, "Contact", pNodeList)
    transpar%ncont = getLength(pNodeList)
    allocate(transpar%contacts(transpar%ncont))
    call readContacts(pNodeList, transpar%contacts, geom, char(buffer), transpar%contactLayerTol,&
        & errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! check for atoms in multiple contact ranges/device or atoms missing from any of these regions
    allocate(atomInRegion(geom%nAtom), source=.false.)
    atomInRegion(transpar%idxdevice(1):transpar%idxdevice(2)) = .true.
    do ii = 1, transpar%nCont
      if (any(atomInRegion(transpar%contacts(ii)%idxrange(1):transpar%contacts(ii)%idxrange(2))))&
          & then
        write(strTmp, "(A,A,A,I0)")"Contact '", trim(transpar%contacts(ii)%name),&
            & "' contains an atom already in the device region or another contact: Atom nr. ",&
            & findloc(atomInRegion(transpar%contacts(ii)%idxrange(1):&
            & transpar%contacts(ii)%idxrange(2)), .true.) + transpar%contacts(ii)%idxrange(1)
        call getItem1(pNodeList, ii, pTmp)
        call detailedError(pTmp, strTmp, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      atomInRegion(transpar%contacts(ii)%idxrange(1):transpar%contacts(ii)%idxrange(2)) = .true.
    end do
    if (any(.not.atomInRegion)) then
      write(strTmp, "(A,I0,A)")"Atom ", findloc(atomInRegion, .false.),&
          & " is not in the device region or any contact"
      call detailedError(root, strTmp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call destroyNodeList(pNodeList)

    transpar%taskUpload = .false.

    select case (char(buffer))
    case ("contacthamiltonian")

      call getChildValue(pTaskType, "ContactId", buffer, errStatus, child=pTmp)
      @:PROPAGATE_ERROR(errStatus)
      call getContactByName(transpar%contacts(:)%name, tolower(trim(unquote(char(buffer)))), pTmp,&
          & contact, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      transpar%taskContInd = contact
      if (.not. geom%tPeriodic) then
        call getChildValue(pTaskType, "ContactSeparation", lateralContactSeparation, errStatus,&
            & 1000.0_dp, modifier=modifier, child=field)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(modifier),lengthUnits,field,lateralContactSeparation, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

      call reduceGeometry(transpar%contacts(contact)%lattice, transpar%contacts(contact)%idxrange,&
          & lateralContactSeparation, geom, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      transpar%ncont = 0

      call getChildValue(root, "writeBinaryContact", transpar%tWriteBinShift, errStatus, .true.)
      @:PROPAGATE_ERROR(errStatus)

    case ("uploadcontacts")

      transpar%taskUpload = .true.

      call getChildValue(root, "readBinaryContact", transpar%tReadBinShift, errStatus, .true.)
      @:PROPAGATE_ERROR(errStatus)

    case default

      call getNodeHSDName(pTaskType, buffer)
      call detailedError(pTask, "Invalid task '" // char(buffer) // "'", errStatus)
      @:PROPAGATE_ERROR(errStatus)

   end select

   call destroyNodeList(pNodeList)

  end subroutine readTransportGeometry


  !> Reduce the geometry for the contact calculation
  subroutine reduceGeometry(contactVec, contactRange, lateralContactSeparation, geom, errStatus)

    !> Vector between principle layers in the contact
    real(dp), intent(in) :: contactVec(3)

    !> Range of atoms in the contact
    integer, intent(in) :: contactRange(2)

    !> Lateral separation distance between contacts in a periodic box
    real(dp), intent(in) :: lateralContactSeparation

    !> atomic geometry
    type(TGeometry), intent(inout) :: geom

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    real(dp) :: contUnitVec(3), dots(3), newLatVecs(3, 3), newOrigin(3)
    real(dp) :: minProj, maxProj
    logical :: mask(3)
    integer :: ind, indPrev, indNext, ii

    if (geom%tPeriodic) then
      contUnitVec = contactVec / sqrt(sum(contactVec**2, dim=1))
      dots = abs(matmul(contUnitVec, geom%latVecs))
      mask = (abs(dots - sqrt(sum(geom%latVecs, dim=1)**2)) < 1e-8_dp)
      if (count(mask) /= 1) then
        @:RAISE_ERROR(errStatus, -1, "Too many lattice vectors parallel to the contact")
      end if
      ! Workaround for bug in Intel compiler (can not use index function)
      ind = 1
      do while (.not. mask(ind))
        ind = ind + 1
      end do
      newLatVecs = geom%latVecs
      newLatVecs(:,ind) = 2.0_dp * contactVec
      newOrigin = geom%origin
    else
      newLatVecs(:,1) = 2.0_dp * contactVec
      mask = abs(contactVec) > 1e-8_dp
      ! Workaround for bug in Intel compiler (can not use index function)
      ind = 1
      do while (.not. mask(ind))
        ind = ind + 1
      end do
      ! Note: ind is one-based, subtract 1 before modulo and add 1 after.
      indNext = modulo(ind + 1 - 1, 3) + 1
      indPrev = modulo(ind - 1 - 1, 3) + 1
      newLatVecs(indNext, 2) = -newLatVecs(ind, 1)
      newLatVecs(ind, 2) = newLatVecs(indNext, 1)
      newLatVecs(indPrev, 2) = 0.0_dp
      newLatVecs(:,3) = cross3(newLatVecs(:,1), newLatVecs(:,2))
      newLatVecs(:,2) = newLatVecs(:,2) / sqrt(sum(newLatVecs(:,2)**2))
      newLatVecs(:,3) = newLatVecs(:,3) / sqrt(sum(newLatVecs(:,3)**2))
      newOrigin(:) = 0.0_dp
    end if
    call reduce(geom, contactRange(1), contactRange(2))
    if (.not. geom%tPeriodic) then
      do ii = 2, 3
        minProj = minval(matmul(newLatVecs(:,ii), geom%coords))
        maxProj = maxval(matmul(newLatVecs(:,ii), geom%coords))
        newLatVecs(:,ii) = ((maxProj - minProj) + lateralContactSeparation) * newLatVecs(:,ii)
      end do
    end if
    call setLattice(geom, newOrigin, newLatVecs)

  end subroutine reduceGeometry


  !> Reads settings for the first layer atoms in principal layers
  subroutine readFirstLayerAtoms(pnode, pls, npl, idxdevice, errStatus, check)

    type(fnode), pointer, intent(in) :: pnode

    !> Start atoms in the principal layers
    integer, allocatable, intent(out) :: pls(:)

    !> Number of principal layers
    integer, intent(out) :: npl

    !> Atoms range of the device
    integer, intent(in) :: idxdevice(2)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> Optional setting to turn on/off check (defaults to on if absent)
    logical, optional, intent(in) :: check


    type(TListInt) :: li
    logical :: checkidx

    checkidx = .true.
    if (present(check)) checkidx = check

    if (associated(pnode)) then
        call init(li)
        call getChildValue(pnode, "", li, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        npl = len(li)
        allocate(pls(npl))
        call asArray(li, pls)
        call destruct(li)
        if (checkidx) then
          if (any(pls < idxdevice(1) .or. &
                  pls > idxdevice(2))) then
            call detailedError(pnode, "First layer atoms must be between "// i2c(idxdevice(1))&
                & // " and " // i2c(idxdevice(2)) // ".", errStatus)
             @:PROPAGATE_ERROR(errStatus)
          end if
        end if
      else
         npl = 1
         allocate(pls(npl))
         pls = (/ 1 /)
      end if

  end subroutine readFirstLayerAtoms


  !> Reads Green's function settings
  subroutine readGreensFunction(pNode, greendens, transpar, tempElec, errStatus)

    !> Input tree
    type(fnode), pointer :: pTmp

    !> Settings for Green's function solver
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

    !> Transport solver settings
    type(TTransPar), intent(inout) :: transpar

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: pNode
    type(fnode), pointer :: field, child1, child2
    real(dp) :: Estep
    integer :: defValue, ii
    type(string) :: buffer, modifier
    logical :: realAxisConv, equilibrium

    type(TListReal) :: fermiBuffer

    greendens%defined = .true.

    if (.not. transpar%defined) then
      !! Fermi level: in case of collinear spin we accept two values
      !! (up and down)
      call init(fermiBuffer)
      call getChildValue(pNode, "FermiLevel", fermiBuffer, errStatus, modifier=modifier)
      @:PROPAGATE_ERROR(errStatus)
      if (len(fermiBuffer) .eq. 1) then
        call asArray(fermiBuffer, greendens%oneFermi)
        greendens%oneFermi(2) = greendens%oneFermi(1)
      else if ( len(fermiBuffer) .eq. 2) then
        call asArray(fermiBuffer, greendens%oneFermi)
      else
        call detailedError(pNode, "FermiLevel accepts 1 or 2 (for collinear spin) values",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call destruct(fermiBuffer)
      call convertUnitHsd(char(modifier), energyUnits, pNode, greendens%oneFermi, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChild(pNode, "FirstLayerAtoms", pTmp, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      call readFirstLayerAtoms(pTmp, greendens%PL, greendens%nPLs, transpar%idxdevice, errStatus,&
          & check=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (.not.associated(pTmp)) then
        call setChildValue(pNode, "FirstLayerAtoms", greendens%PL, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      !call getChild(pNode, "ContactPLs", pTmp, errStatus, requested=.false.)
      !if (associated(pTmp)) then
      !  call init(li)
      !  call getChildValue(pTmp, "", li, errStatus)
      !  @:PROPAGATE_ERROR(errStatus)
      !  allocate(transpar%cblk(len(li)))
      !  call asArray(li,transpar%cblk)
      !  call destruct(li)
      !end if
      allocate(greendens%kbT(1))
      greendens%kbT(:) = tempElec
    else
      if (transpar%ncont > 0) then
        allocate(greendens%kbT(transpar%ncont))
        do ii = 1, transpar%ncont
          if (transpar%contacts(ii)%kbT .ge. 0.0_dp) then
            greendens%kbT(ii) = transpar%contacts(ii)%kbT
          else
            greendens%kbT(ii) = tempElec
          end if
        enddo
      end if
    end if

    call getChildValue(pNode, "LocalCurrents", greendens%doLocalCurr, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "Verbosity", greendens%verbose, errStatus, 51)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "Delta", greendens%delta, errStatus, 1.0e-5_dp, modifier=modifier,&
        & child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), energyUnits, field, greendens%delta, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "ReadSurfaceGFs", greendens%readSGF, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "SaveSurfaceGFs", greendens%saveSGF, errStatus,&
        & .not.greendens%readSGF)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "ContourPoints", greendens%nP(1:2), errStatus, [ 20, 20 ])
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "EnclosedPoles",  greendens%nPoles, errStatus, 3)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "LowestEnergy", greendens%enLow, errStatus, -2.0_dp,&
        & modifier=modifier, child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), energyUnits, field, greendens%enLow, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "FermiCutoff", greendens%nkT, errStatus, 10)
    @:PROPAGATE_ERROR(errStatus)
    ! Fermi energy had not been set by other means yet

    ! Non equilibrium integration along real axis:
    ! The code will perform the integration if the number of points is larger
    ! than zero, no matter if there's bias or not.
    ! Therefore I restored the default on the energy step, as it works at zero
    ! bias and it scales flawlessy with increasing bias
    ! It is still allowed to directly set the number of points, if preferred
    ! libNEGF only wants the number of points in input
    call getChild(pNode, "RealAxisPoints", child1, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChild(pNode, "RealAxisStep", child2, errStatus, requested=.false., modifier=buffer)
    @:PROPAGATE_ERROR(errStatus)
      realAxisConv = .false.
      ! Set a bool to verify if all contacts are at the same potential (if so,
      ! no points are needed)
      equilibrium = .true.
      do ii = 2, transpar%ncont
        if (transpar%contacts(1)%potential .ne. transpar%contacts(ii)%potential &
           & .or. transpar%contacts(1)%kbT .ne. transpar%contacts(ii)%kbT ) then
           equilibrium = .false.
        end if
      end do

      ! Both Points and Step cannot be specified
      if  (associated (child1) .and. associated(child2)) then
        call detailedError(child1, "RealAxisPoints and RealAxisStep cannot be specified together.",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      ! If only one is specified, take it as valid value
      else if (associated(child1)) then
        call getChildValue(pNode, "RealAxisPoints", greendens%nP(3), errStatus)
        @:PROPAGATE_ERROR(errStatus)
      else if (associated(child2)) then
        call getChildValue(pNode, "RealAxisStep", Estep, errStatus, child=child2,&
            & modifier=modifier)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(modifier), energyUnits, child2, Estep, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        realAxisConv = .true.
      ! If the system is under equilibrium we set the number of
      ! points to zero
      else if (equilibrium) then
        call getChildValue(pNode, "RealAxisPoints", greendens%nP(3), errStatus, 0, child=child1)
        @:PROPAGATE_ERROR(errStatus)
      else
        !Default is a point every 1500H
        call getChildValue(pNode, "RealAxisStep", Estep, errStatus, 6.65e-4_dp, modifier=modifier,&
            & child=child2)
        @:PROPAGATE_ERROR(errStatus)
        realAxisConv = .true.
      end if
      ! RealAxisConv means that we have a step and we convert it in a number
      ! of points
      if (realAxisConv) then
        defValue = int(1.0_dp/Estep &
          & * (maxval(transpar%contacts(:)%potential) &
          & - minval(transpar%contacts(:)%potential) + &
          & 2 * greendens%nKT * maxval(greendens%kbT)))
        greendens%nP(3) = defvalue
        ! call getChildValue(pNode, "RealAxisPoints", greendens%nP(3), errStatus, defvalue,&
        !     & child=child1)
        ! @:PROPAGATE_ERROR(errStatus)
      end if

  end subroutine readGreensFunction
#:endif


#:if WITH_POISSON

  !> Read in Poisson related data
#:if WITH_TRANSPORT
  subroutine readPoisson(pNode, poisson, tPeriodic, transpar, latVecs, updateSccAfterDiag,&
      & errStatus)
#:else
  subroutine readPoisson(pNode, poisson, tPeriodic, latVecs, updateSccAfterDiag, errStatus)
#:endif

    !> Input tree
    type(fnode), pointer :: pNode

    !> data type for Poisson solver settings
    type(TPoissonInfo), intent(inout) :: poisson

    !> Is this a periodic calculation
    logical, intent(in) :: tPeriodic

  #:if WITH_TRANSPORT
    !> Parameters of the transport calculation
    type(TTransPar), intent(inout) :: transpar
  #:endif

    !> Lattice vectors if periodic
    real(dp), allocatable, intent(in) :: latVecs(:,:)

    !> Whether Scc should be updated with the output charges (obtained after diagonalisation)
    logical, intent(out) :: updateSccAfterDiag

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: pTmp, pTmp2, pChild, field
    type(string) :: buffer, modifier
    real(dp) :: denstol, gatelength_l
    logical :: needsPoissonBox

  #:if WITH_TRANSPORT
    integer :: ii
  #:endif

    poisson%defined = .true.
    needsPoissonBox = .not. tPeriodic
  #:if WITH_TRANSPORT
    needsPoissonBox = needsPoissonBox .or. transpar%tPeriodic1D .or. transpar%nCont == 1
  #:endif

    if (needsPoissonBox) then
    #:if WITH_TRANSPORT
      if (transpar%nCont == 1 .and. .not. transpar%tPeriodic1D) then
        poisson%poissBox(:) = 0.0_dp
        do ii = 1, 3
          if (ii == transpar%contacts(1)%dir) then
            call getChildValue(pNode, "PoissonThickness", poisson%poissBox(ii), errStatus,&
                & modifier=modifier, child=field)
            @:PROPAGATE_ERROR(errStatus)
            call convertUnitHsd(char(modifier), lengthUnits, field, poisson%poissBox, errStatus)
            @:PROPAGATE_ERROR(errStatus)
          else
            poisson%poissBox(ii) = sqrt(sum(latVecs(:,ii)**2))
          end if
        end do
      else
        call getChildValue(pNode, "PoissonBox", poisson%poissBox, errStatus, modifier=modifier,&
            & child=field)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(modifier), lengthUnits, field, poisson%poissBox, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    #:else
      call getChildValue(pNode, "PoissonBox", poisson%poissBox, errStatus, modifier=modifier,&
          & child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%poissBox, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:endif
    end if

    poisson%foundBox = needsPoissonBox
    call getChildValue(pNode, "MinimalGrid", poisson%poissGrid, errStatus,&
        & [0.3_dp, 0.3_dp, 0.3_dp], modifier=modifier, child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), lengthUnits, field, poisson%poissGrid, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "NumericalNorm", poisson%numericNorm, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChild(pNode, "AtomDensityCutoff", pTmp, errStatus, requested=.false., modifier=modifier)
    @:PROPAGATE_ERROR(errStatus)
    call getChild(pNode, "AtomDensityTolerance", pTmp2, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(pTmp) .and. associated(pTmp2)) then
      call detailedError(pNode, "Either one of the tags AtomDensityCutoff or AtomDensityTolerance&
          & can be specified.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else if (associated(pTmp)) then
      call getChildValue(pTmp, "", poisson%maxRadAtomDens, errStatus, default=14.0_dp,&
          & modifier=modifier)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, pTmp, poisson%maxRadAtomDens, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (poisson%maxRadAtomDens <= 0.0_dp) then
        call detailedError(pTmp2, "Atom density cutoff must be > 0", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    else
      call getChildValue(pNode, "AtomDensityTolerance", denstol, errStatus, 1e-6_dp, child=pTmp2)
      @:PROPAGATE_ERROR(errStatus)
      if (denstol <= 0.0_dp) then
        call detailedError(pTmp2, "Atom density tolerance must be > 0", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      ! Negative value to signal automatic determination
      poisson%maxRadAtomDens = -denstol
    end if

    call getChildValue(pNode, "CutoffCheck", poisson%cutoffcheck, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "Verbosity", poisson%verbose, errStatus, 51)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "SavePotential", poisson%savePotential, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "PoissonAccuracy", poisson%poissAcc, errStatus, 1.0e-6_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "BuildBulkPotential", poisson%bulkBC, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "ReadOldBulkPotential", poisson%readBulkPot, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "RecomputeAfterDensity", updateSccAfterDiag, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(pNode, "MaxPoissonIterations", poisson%maxPoissIter, errStatus, 60)
    @:PROPAGATE_ERROR(errStatus)

    poisson%overrideBC(:) = 0
    call getChild(pNode, "OverrideDefaultBC", pTmp, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(pTmp)) then
      call getPoissonBoundaryConditionOverrides(pTmp, [ 1, 2 ], poisson%overrideBC, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getChildValue(pNode, "OverrideBulkBC", pTmp, errStatus, "none")
    @:PROPAGATE_ERROR(errStatus)
    poisson%overrBulkBC(:) = -1
    if (associated(pNode)) then
      call getPoissonBoundaryConditionOverrides(pTmp, [ 0, 1, 2 ], poisson%overrBulkBC, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getChildValue(pNode, "BoundaryRegion", pTmp, errStatus, "global")
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(pTmp, buffer)
    select case(char(buffer))
    case ("global")
      poisson%localBCType = "G"
    case ("square")
      poisson%localBCType = "S"
      call getChildValue(pTmp, "BufferLength", poisson%bufferLocBC, errStatus, 9.0_dp,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%bufferLocBC, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("circle")
      poisson%localBCType = "C"
      call getChildValue(pTmp, "BufferLength", poisson%bufferLocBC, errStatus, 9.0_dp,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%bufferLocBC, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case default
      call getNodeHSDName(pTmp, buffer)
      call detailedError(pTmp, "Invalid boundary region type '" // char(buffer) // "'", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    call getChildValue(pNode, "BoxExtension", poisson%bufferBox, errStatus, 0.0_dp,&
        & modifier=modifier, child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), lengthUnits, field, poisson%bufferBox, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (poisson%bufferBox.lt.0.0_dp) then
      call detailedError(pNode, "BoxExtension must be a positive number", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    endif

    ! PARSE GATE OPTIONS
    call getChildValue(pNode, "Gate", pTmp2, errStatus, "none", child=pChild)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(pTmp2, buffer)

    poisson%insLength = 0.0_dp
    poisson%insRad = 0.0_dp
    select case(char(buffer))
    case ("none")
      poisson%gateType = "N"
    case ("planar")
      poisson%gateType = "P"
      call getChildValue(pTmp2, "GateLength", poisson%gateLength_l, errStatus, 0.0_dp,&
          & modifier= modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateLength_l, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      gatelength_l = poisson%gateLength_l !avoids a warning on intents
      call getChildValue(pTmp2, "GateLength_l", poisson%gateLength_l, errStatus, gateLength_l,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateLength_l, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(pTmp2, "GateLength_t", poisson%gateLength_t, errStatus,&
          & poisson%gateLength_l, modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateLength_t, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(pTmp2, "GateDistance", poisson%gateRad, errStatus, 0.0_dp,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateRad, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(pTmp2, "GatePotential", poisson%gatepot, errStatus, 0.0_dp,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, field, poisson%gatepot, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      ! call getChildValue(pTmp2, "GateDirection", poisson%gatedir, errStatus, 2)
      ! @:PROPAGATE_ERROR(errStatus)
      poisson%gatedir = 2

    case ("cylindrical")
      poisson%gateType = "C"
      call getChildValue(pTmp2, "GateLength", poisson%gateLength_l, errStatus, 0.0_dp,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateLength_l, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(pTmp2, "GateRadius", poisson%gateRad, errStatus, 0.0_dp,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateRad, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(pTmp2, "GatePotential", poisson%gatepot, errStatus, 0.0_dp,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, field, poisson%gatepot, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case default
      call getNodeHSDName(pTmp2, buffer)
      call detailedError(pTmp2, "Invalid gate type '" // char(buffer) // "'", errStatus)
      @:PROPAGATE_ERROR(errStatus)

    end select

    call getChildValue(pNode, "MaxParallelNodes", poisson%maxNumNodes, errStatus, 1)
    @:PROPAGATE_ERROR(errStatus)

    poisson%scratch = "contacts"

  end subroutine readPoisson


  !> Over-rides the boundary conditions on the Poisson solver
  subroutine getPoissonBoundaryConditionOverrides(pNode, availableConditions, overrideBC, errStatus)

    !> Input data tree
    type(fnode), pointer, intent(in) :: pNode

    !> List of conditions that can be set as choices
    integer, intent(in) :: availableConditions(:)

    !> Array of boundary condition types on the 6 faces of the box, 0 for use of default
    integer, intent(inout) :: overrideBC(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer, parameter :: PERIODIC_BC = 0
    integer, parameter :: DIRICHLET_BC = 1
    integer, parameter :: NEUMANN_BC = 2
    character(10), parameter :: bcstr(0:2) = &
        & [ character(10) :: "Periodic", "Dirichlet", "Neumann" ]
    integer :: bctype, iBC
    integer :: faceBC, oppositeBC
    integer :: ii
    type(TListString) :: lStr
    type(fnode), pointer :: pNode2, pChild
    character(lc) :: strTmp

    do iBC = 1, size(availableConditions)
      bctype = availableConditions(iBC)
      call getChild(pNode, trim(bcstr(bctype)), pNode2, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(pNode2)) then
        call init(lStr)
        call getChildValue(pNode2, "boundaries", lStr, errStatus, child=pChild)
        @:PROPAGATE_ERROR(errStatus)
        if (len(lStr).gt.6) then
          call detailedError(pChild,"boundaries must be 6 or less", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        do ii = 1, len(lStr)
          call get(lStr, strTmp, ii)
          select case(trim(strTmp))
          case("x")
            overrideBC(1) = bctype
            overrideBC(2) = bctype
          case("xmin")
            overrideBC(1) = bctype
          case("xmax")
            overrideBC(2) = bctype
          case("y")
            overrideBC(3) = bctype
            overrideBC(4) = bctype
          case("ymin")
            overrideBC(3) = bctype
          case("ymax")
            overrideBC(4) = bctype
          case("z")
            overrideBC(5) = bctype
            overrideBC(6) = bctype
          case("zmin")
            overrideBC(5) = bctype
          case("zmax")
            overrideBC(6) = bctype
          end select
        end do
        call destruct(lStr)
      end if
    end do

    ! If face is set to periodic, opposite one should be the same
    do ii = 1, 3
      faceBC = overrideBC(2 * ii)
      oppositeBC = overrideBC(2 * ii - 1)
      if (faceBC == PERIODIC_BC &
          & .and. oppositeBC /= faceBC) then
        call detailedError(pChild, "Periodic override must be set both min max", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end do

  end subroutine getPoissonBoundaryConditionOverrides

#:endif


#:if WITH_TRANSPORT
  !> Sanity checking of atom ranges and returning contact vector and direction.
  subroutine getContactVector(atomrange, geom, id, name, pContact, contactLayerTol, contactVec,&
      & contactDir, errStatus)

    !> Range of atoms in the contact
    integer, intent(in) :: atomrange(2)

    !> Atomic geometry, including the contact atoms
    type(TGeometry), intent(in) :: geom

    !> Index for this contact
    integer, intent(in) :: id

    !> Contact name
    character(mc), intent(in) :: name

    !> Node in the parser, needed for error handling
    type(fnode), pointer :: pContact

    !> Allowed discrepancy in positions of atoms between the contact's two  principle layers
    real(dp), intent(in) :: contactLayerTol

    !> Vector direction between principal layers in the contact
    real(dp), intent(out) :: contactVec(3)

    !> Which supercell vector the contact vector is parallel to
    integer, intent(out) :: contactDir

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: iStart, iStart2, iEnd, ii
    logical :: mask(3)
    character(lc) :: errorStr

    !! Sanity check for the atom ranges
    iStart = atomrange(1)
    iEnd = atomrange(2)
    if (iStart < 1 .or. iEnd < 1 .or. iStart > geom%nAtom .or. iEnd > geom%nAtom) then
      call detailedError(pContact, "Invalid atom range '" // i2c(iStart) &
          &// " " // i2c(iEnd) // "', values should be between " // i2c(1) &
          &// " and " // i2c(geom%nAtom) // ".", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    if (iEnd < iStart) then
      call detailedError(pContact, "Invalid atom order in contact '" // i2c(iStart) // " " //&
          & i2c(iEnd) // "', should be asscending order.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    if (mod(iEnd - iStart + 1, 2) /= 0) then
      call detailedError(pContact, "Nr. of atoms in the contact must be even", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! Determining intra-contact layer vector
    iStart2 = iStart + (iEnd - iStart + 1) / 2
    contactVec = geom%coords(:,iStart) - geom%coords(:,iStart2)

    if (any(sum( (geom%coords(:,iStart:iStart2-1) - geom%coords(:,iStart2:iEnd)&
        & - spread(contactVec, dim=2, ncopies=iStart2-iStart))**2, dim=1) > contactLayerTol**2))&
        & then
      write(stdout,"(1X,A,I0,A,I0)")'Contact vector defined from atoms ', iStart, ' and ',iStart2
      write(stdout,"(1X,A,I0,'-',I0)")'Contact layer 1 atoms: ',iStart, iStart2-1
      write(stdout,"(1X,A,I0,'-',I0)")'Contact layer 2 atoms: ',iStart2, iEnd
      do ii = 0, iStart2 -1 -iStart
        if (sum((geom%coords(:,ii+iStart)-geom%coords(:,ii+iStart2) - contactVec)**2)&
            & > contactLayerTol**2) then
          write(stdout,"(1X,A,I0,A,I0,A)")'Atoms ',iStart+ii, ' and ', iStart2+ii,&
              & ' inconsistent with the contact vector.'
          exit
        end if
      end do
      write(stdout,*)'Mismatches in atomic positions in the two layers:'
      write(stdout,"(3F20.12)")((geom%coords(:,iStart:iStart2-1) - geom%coords(:,iStart2:iEnd)&
          & - spread(contactVec(:), dim=2, ncopies=iStart2-iStart))) * Bohr__AA

      write (errorStr,"('Contact ',A,' (',A,') does not consist of two rigidly shifted layers')")&
          & i2c(id), trim(name)
      @:RAISE_ERROR(errStatus, -1, errorStr)

    end if

    ! Determine to which axis the contact vector is parallel.
    mask = (abs(abs(contactVec) - sqrt(sum(contactVec**2))) < 1.0e-8_dp)
    if (count(mask) /= 1) then
      call warning("Contact vector " // i2c(id) // " not parallel to any of the coordinate axis.")
      contactDir = 0
    else
      ! Workaround for bug in Intel compiler (can not use index function)
      contactDir = 1
      do while (.not. mask(contactDir))
        contactDir = contactDir + 1
      end do
    end if

  end subroutine getContactVector


  !> Read dephasing block
  subroutine readDephasing(node, orb, geom, tp, tundos, errStatus)

    !> Input tree node
    type(fnode), pointer :: node

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Atomic geometry, including the contact atoms
    type(TGeometry), intent(in) :: geom

    !> Parameters of the transport calculation
    type(TTransPar), intent(inout) :: tp

    !> Parameters of tunneling and dos calculation
    type(TNEGFTunDos), intent(inout) :: tundos

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child

    call getChild(node, "VibronicElastic", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child)) then
      tp%tDephasingVE = .true.
      call readElPh(child, tundos%elph, geom, orb, tp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getChildValue(node, "BuettikerProbes", value1, errStatus, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(value1)) then
      tp%tDephasingBP = .true.
      call readDephasingBP(child, tundos%bp, geom, orb, tp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    ! Lowdin transformations involve dense matrices and works only in small systems
    ! For the dftb+ official release the options are disabled
    tp%tOrthonormal = .false.
    tp%tOrthonormalDevice = .false.
    ! call getChildValue(node, "Orthonormal", tp%tOrthonormal, errStatus, .false.)
    ! @:PROPAGATE_ERROR(errStatus)
    ! call getChildValue(node, "OrthonormalDevice", tp%tOrthonormalDevice, errStatus, .false.)
    ! @:PROPAGATE_ERROR(errStatus)
    tp%tNoGeometry = .false.
    tp%NumStates = 0

  end subroutine readDephasing


  !> Read Electron-Phonon blocks (for density and/or current calculation)
  subroutine readElPh(node, elph, geom, orb, tp, errStatus)

    !> Input node in the tree
    type(fnode), pointer :: node

    !> container for electron-phonon parameters
    type(TElPh), intent(inout) :: elph

    !> Geometry type
    type(TGeometry), intent(in) :: geom

    !> Orbitals infos
    type(TOrbitals), intent(in) :: orb

    !> Transport parameter type
    type(TTransPar), intent(in) :: tp

    !> Error status
    type(TStatus), intent(inout) :: errStatus


    logical :: block_model, semilocal_model

    elph%defined = .true.
    !! Only local el-ph model is defined (elastic for now)
    elph%model = 1

    call getChildValue(node, "MaxSCBAIterations", elph%scba_niter, errStatus, default=100)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "atomBlock", block_model, errStatus, default=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (block_model) then
      elph%model = 2
    endif

    !BUG: semilocal model crashes because of access of S before its allocation
    !     this because initDephasing was moved into initprogram
    call getChildValue(node, "semiLocal", semilocal_model, errStatus, default=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (semilocal_model) then
      call detailedError(node, "Semilocal dephasing causes crash and has been "//&
          & "temporarily disabled", errStatus)
      @:PROPAGATE_ERROR(errStatus)
      elph%model = 3
    endif

    call readCoupling(node, elph, geom, orb, tp, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readElPh


  !> Read Buettiker probe dephasing blocks (for density and/or current calculation)
  subroutine readDephasingBP(node, elph, geom, orb, tp, errStatus)

    !> Node in input document tree
    type(fnode), pointer :: node

    !> container for buttiker-probes parameters
    type(TElPh), intent(inout) :: elph

    !> Geometry type
    type(TGeometry), intent(in) :: geom

    !> Orbitals infos
    type(TOrbitals), intent(in) :: orb

    !> Transport parameter type
    type(TTransPar), intent(inout) :: tp

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    logical :: block_model, semilocal_model
    type(string) :: model
    type(fnode), pointer :: dephModel

    call detailedError(node, "Buettiker probes are still under development", errStatus)
    @:PROPAGATE_ERROR(errStatus)

    elph%defined = .true.
    call getChildValue(node, "", dephModel, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName2(dephModel, model)

    select case(char(model))
    case("dephasingprobes")
      !! Currently only zeroCurrent condition is implemented
      !! This corresponds to elastic dephasing probes
      tp%tZeroCurrent=.true.
      !! Only local bp model is defined (elastic for now)
    case("voltageprobes")
      call detailedError(dephModel,"voltageProbes have been not implemented yet", errStatus)
      @:PROPAGATE_ERROR(errStatus)
      tp%tZeroCurrent=.false.
    case default
      call detailedError(dephModel, "unknown model", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    elph%model = 1

    call getChildValue(dephModel, "MaxSCBAIterations", elph%scba_niter, errStatus, default=100)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(dephModel, "atomBlock", block_model, errStatus, default=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (block_model) then
      elph%model = 2
    endif

    !BUG: semilocal model crashes because of access of S before its allocation
    !     this because initDephasing occurs in initprogram
    call getChildValue(dephModel, "semiLocal", semilocal_model, errStatus, default=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (semilocal_model) then
      call detailedError(dephModel, "semilocal dephasing is not working yet", errStatus)
      @:PROPAGATE_ERROR(errStatus)
      elph%model = 3
    endif

    call readCoupling(dephModel, elph, geom, orb, tp, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readDephasingBP


  !> Reads coupling strength and mode for dephasing
  !> 2 modes support, constant or specified per each orbital
  subroutine readCoupling(node, elph, geom, orb, tp, errStatus)

    !> Node in the input tree
    type(fnode), pointer :: node

    !> container for buttiker-probes parameters
    type(TElPh), intent(inout) :: elph

    !> Geometry type
    type(TGeometry), intent(in) :: geom

    !> Orbitals infos
    type(TOrbitals), intent(in) :: orb

    !> Transport parameter type
    type(TTransPar), intent(in) :: tp

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(string) :: buffer, method, modifier, modifier2
    type(fnode), pointer :: val, child, child2, child3, child4, field
    type(fnodeList), pointer :: children
    integer :: norbs, ii, jj, iAt
    integer :: atm_range(2)
    real(dp) :: rTmp
    integer, allocatable :: tmpI1(:)
    real(dp), allocatable :: atmCoupling(:)

    !! Allocate coupling array
    norbs = 0
    if (tp%defined) then
      atm_range(1) = tp%idxdevice(1)
      atm_range(2) = tp%idxdevice(2)
    else
      atm_range(1) = 1
      atm_range(2) = geom%nAtom
    endif
    do ii=atm_range(1), atm_range(2)
      norbs = norbs + orb%nOrbAtom(ii)
    enddo
    allocate(elph%coupling(norbs))
    elph%coupling(:) = 0.d0

    elph%orbsperatm = orb%nOrbAtom(atm_range(1):atm_range(2))

    call getChildValue(node, "Coupling", val, errStatus, "", child=child,&
        & allowEmptyValue=.true., modifier=modifier, dummyValue=.true., list=.false.)
    @:PROPAGATE_ERROR(errStatus)

    call getNodeName(val, method)

    ! This reads also things like:  "Coupling [eV] = 0.34"
    !if (is_numeric(char(method))) then
    !  call getChildValue(node, "Coupling", rTmp, errStatus, child=field)
    !  @:PROPAGATE_ERROR(errStatus)
    !  call convertUnitHsd(char(modifier), energyUnits, field, rTmp, errStatus)
    !  elph%coupling = rTmp
    !  return
    !end if

    select case (char(method))
    case ("allorbitals")
      call getChild(child, "AllOrbitals", child2, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(child2, "", elph%coupling, errStatus, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, field, elph%coupling, errStatus)
      @:PROPAGATE_ERROR(errStatus)

    case ("atomcoupling")
      call getChild(child, "AtomCoupling", child2, errStatus, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      allocate(atmCoupling(atm_range(2)-atm_range(1)+1))
      atmCoupling = 0.d0
      call getChildren(child2, "AtomList", children)
      do ii = 1, getLength(children)
        call getItem1(children, ii, child3)
        call getChildValue(child3, "Atoms", buffer, errStatus, child=child4, multiple=.true.)
        @:PROPAGATE_ERROR(errStatus)
        call getSelectedAtomIndices(child4, char(buffer), geom%speciesNames, geom%species, tmpI1,&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call getChildValue(child3, "Value", rTmp, errStatus, child=field, modifier=modifier2)
        @:PROPAGATE_ERROR(errStatus)
        ! If not defined, use common unit modifier defined after Coupling
        if (len(modifier2)==0) then
          call convertUnitHsd(char(modifier), energyUnits, field, rTmp, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        else
          call convertUnitHsd(char(modifier2), energyUnits, field, rTmp, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
        do jj=1, size(tmpI1)
          iAt = tmpI1(jj)
          if (atmCoupling(iAt) /= 0.0_dp) then
            call detailedWarning(child3, "Previous setting of coupling &
                &for atom" // i2c(iAt) // " has been overwritten")
          end if
          atmCoupling(iAt) = rTmp
        end do
      end do
      call destroyNodeList(children)

      ! Transform atom coupling in orbital coupling
      norbs = 0
      do ii=atm_range(1), atm_range(2)
        elph%coupling(norbs + 1:norbs + orb%nOrbAtom(ii)) = atmCoupling(ii)
        norbs = norbs + orb%nOrbAtom(ii)
      enddo
      deallocate(atmCoupling)

    case ("constant")
      call getChildValue(child, "Constant", rtmp, errStatus, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, field, rTmp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      elph%coupling = rTmp

    case default
      call detailedError(node, "Coupling definition unknown", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readCoupling


  !> Read Tunneling and Dos options from analysis block
  subroutine readTunAndDos(root, orb, geo, tundos, transpar, tempElec, errStatus)
    type(fnode), pointer :: root
    type(TOrbitals), intent(in) :: orb
    type(TGeometry), intent(in) :: geo

    !> tundos is the container to be filled
    type(TNEGFTunDos), intent(inout) :: tundos
    type(TTransPar), intent(inout) :: transpar
    real(dp), intent(in) :: tempElec

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: pTmp, pNode, field
    type(fnodeList), pointer :: pNodeList
    integer :: ii, jj, ind, ncont, nKT
    real(dp) :: eRange(2), eRangeDefault(2)
    type(string) :: modifier
    type(TWrappedInt1), allocatable :: iAtInRegion(:)
    logical, allocatable :: tShellResInRegion(:)
    character(lc), allocatable :: regionLabelPrefixes(:)
    type(TListReal) :: temperature

    tundos%defined = .true.

    ! ncont is needed for contact option allocation
    ncont = transpar%ncont

    call getChildValue(root, "Verbosity", tundos%verbose, errStatus, 51)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(root, "WriteLDOS", tundos%writeLDOS, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(root, "WriteTunn", tundos%writeTunn, errStatus, .true.)
    @:PROPAGATE_ERROR(errStatus)

    ! Read Temperature. Can override contact definition
    allocate(tundos%kbT(ncont))
    call getChild(root, "ContactTemperature", pTmp, errStatus, modifier=modifier, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(pTmp)) then
      call init(temperature)
      call getChildValue(pTmp, "", temperature, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (len(temperature) .ne. ncont) then
        call detailedError(root, "ContactTemperature does not match the number of contacts",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call asArray(temperature, tundos%kbT)
      call destruct(temperature)
      call convertUnitHsd(char(modifier), energyUnits, pTmp, tundos%kbT, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else
      do ii = 1, ncont
        if (transpar%contacts(ii)%kbT >= 0) then
          tundos%kbT(ii) = transpar%contacts(ii)%kbT
        else
          tundos%kbT(ii) = tempElec
        end if
      end do
    end if

    ! Parsing of energy range
    ! If the calculation is in equilibrium (all potentials to 0.0)
    ! then an energy range and step must be specified (it is assumed
    ! that the user use this filed to calculate a DOS or T(E) )
    ! If the calculation is out of equilibrium, a default similar to
    ! GreensFunction RealAxisStep is set to ensure that the current
    ! can be calculated without manually specify the energy parameters.

    if (all(transpar%contacts(:)%potential.eq.0.0)) then
      ! No default meaningful
      call getChildValue(root, "EnergyRange", eRange, errStatus, modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, field, eRange, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(root, "EnergyStep", tundos%estep, errStatus, modifier=modifier,&
          & child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, field, tundos%estep, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    else
      ! Default meaningful
      ! nKT is set to GreensFunction default, i.e. 10
      ! I avoid an explicit nKT option because I find it confusing here
      ! (it makes sense only out of equilibrium)
      ! Emin = min(-mu); Emax=max(-mu) where mu is Vi-min(Efi)
      ! Note: if Efi != min(Efi) a built in potential is added in poisson
      ! to aling the leads, we don't need to include it here
      nKT = 10
      eRangeDefault(1) = minval(-1.0*transpar%contacts(:)%potential) + &
                        & minval(1.0*transpar%contacts(:)%eFermi(1)) -   &
                        & nKT * maxval(tundos%kbT)
      eRangeDefault(2) = maxval(-1.0*transpar%contacts(:)%potential) + &
                        & minval(transpar%contacts(:)%eFermi(1)) +   &
                        & nKT * maxval(tundos%kbT)
      call getChildValue(root, "EnergyStep", tundos%estep, errStatus, 6.65e-4_dp,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, field, tundos%estep, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(root, "EnergyRange", eRange, errStatus, eRangeDefault, modifier=modifier,&
          & child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, field, eRange, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    tundos%emin = eRange(1)
    tundos%emax = eRange(2)
    ! Terminal currents
    call getChild(root, "TerminalCurrents", pTmp, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
      if (associated(pTmp)) then
        call getChildren(pTmp, "EmitterCollector", pNodeList)
        allocate(tundos%ni(getLength(pNodeList)))
        allocate(tundos%nf(getLength(pNodeList)))
        do ii = 1, getLength(pNodeList)
          call getItem1(pNodeList, ii, pNode)
          call getEmitterCollectorByName(pNode, tundos%ni(ii), tundos%nf(ii),&
              & transpar%contacts(:)%name, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end do
        call destroyNodeList(pNodeList)
      else
        allocate(tundos%ni(ncont-1) )
        allocate(tundos%nf(ncont-1) )
        call setChild(root, "TerminalCurrents", pTmp, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        ind = 1
        do ii = 1, 1
          do jj = ii + 1, ncont
            call setChildValue(pTmp, "EmitterCollector", [transpar%contacts(ii)%name,&
                & transpar%contacts(jj)%name], errStatus)
            @:PROPAGATE_ERROR(errStatus)
            tundos%ni(ind) = ii
            tundos%nf(ind) = jj
            ind = ind + 1
          end do
        end do
      end if
      call getChildValue(root, "Delta", tundos%delta, errStatus, 1.0e-5_dp, modifier=modifier,&
          & child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, field, tundos%delta, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(root, "BroadeningDelta", tundos%broadeningDelta, errStatus, 0.0_dp,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, field, tundos%broadeningDelta, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call readPDOSRegions(root, geo, transpar%idxdevice, iAtInRegion, tShellResInRegion,&
          & regionLabelPrefixes, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      if (allocated(iAtInRegion)) then
        call transformPdosRegionInfo(iAtInRegion, tShellResInRegion, &
            & regionLabelPrefixes, orb, geo%species, tundos%dosOrbitals, &
            & tundos%dosLabels)
      end if

  end subroutine readTunAndDos


  !> Read bias information, used in Analysis and Green's function eigensolver
  subroutine readContacts(pNodeList, contacts, geom, task, contactLayerTol, errStatus)

    !> Node to process
    type(fnodeList), pointer :: pNodeList

    !> Contacts
    type(ContactInfo), allocatable, dimension(:), intent(inout) :: contacts

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    !> What type of transport-related calculation is this?
    character(*), intent(in) :: task

    !> Tolerance to distortion of contact vectors
    real(dp), intent(out) :: contactLayerTol

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: ii
    type(fnode), pointer :: field, pNode, pTmp, child1, child2
    type(string) :: buffer, modifier
    type(TListReal) :: fermiBuffer

    do ii = 1, size(contacts)

      contacts(ii)%wideBand = .false.
      contacts(ii)%wideBandDos = 0.0_dp

      call getItem1(pNodeList, ii, pNode)
      call getChildValue(pNode, "Id", buffer, errStatus, child=pTmp)
      @:PROPAGATE_ERROR(errStatus)
      buffer = tolower(trim(unquote(char(buffer))))
      if (len(buffer) > mc) then
        call detailedError(pTmp, "Contact id may not be longer than " // i2c(mc) // " characters.",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      contacts(ii)%name = char(buffer)
      if (any(contacts(1:ii-1)%name == contacts(ii)%name)) then
        call detailedError(pTmp, "Contact id '" // trim(contacts(ii)%name) //  "' already in use",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if

      call getChildValue(pNode, "PLShiftTolerance", contactLayerTol, errStatus, 1e-5_dp,&
          & modifier=modifier, child=field)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, field, contactLayerTol, errStatus)
      @:PROPAGATE_ERROR(errStatus)

      call getChildValue(pNode, "AtomRange", contacts(ii)%idxrange, errStatus, child=pTmp)
      @:PROPAGATE_ERROR(errStatus)
      call getContactVector(contacts(ii)%idxrange, geom, ii, contacts(ii)%name, pTmp,&
          & contactLayerTol, contacts(ii)%lattice, contacts(ii)%dir, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      contacts(ii)%length = sqrt(sum(contacts(ii)%lattice**2))

      ! Contact temperatures. A negative default is used so it is quite clear when the user sets a
      ! different value. In such a case this overrides values defined in the Filling block
      call getChild(pNode,"Temperature", field, errStatus, modifier=modifier, requested=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (associated(field)) then
        call getChildValue(pNode, "Temperature", contacts(ii)%kbT, errStatus, 0.0_dp,&
            & modifier=modifier, child=field)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(modifier), energyUnits, field, contacts(ii)%kbT, errStatus)
        @:PROPAGATE_ERROR(errStatus)
      else
        contacts(ii)%kbT = -1.0_dp ! -1.0 simply means 'not defined'
      end if

      if (task .eq. "uploadcontacts") then
        call getChildValue(pNode, "Potential", contacts(ii)%potential, errStatus, 0.0_dp,&
            & modifier=modifier, child=field)
        @:PROPAGATE_ERROR(errStatus)
        call convertUnitHsd(char(modifier), energyUnits, field, contacts(ii)%potential, errStatus)
        @:PROPAGATE_ERROR(errStatus)

        call getChildValue(pNode, "WideBand", contacts(ii)%wideBand, errStatus, .false.)
        @:PROPAGATE_ERROR(errStatus)

        if (contacts(ii)%wideBand) then

          ! WideBandApproximation is defined as energy spacing between levels of the contact. In the
          ! code the inverse value (Density of states) is used. Convert the negf input
          ! value. Default is 20 / e eV.
          call getChildValue(pNode, "LevelSpacing", contacts(ii)%wideBandDos, errStatus, 0.735_dp,&
              & modifier=modifier, child=field)
          @:PROPAGATE_ERROR(errStatus)
          call convertUnitHsd(char(modifier), energyUnits, field, contacts(ii)%wideBandDos,&
              & errStatus)
          @:PROPAGATE_ERROR(errStatus)
          contacts(ii)%wideBandDos = 1.d0 / contacts(ii)%wideBandDos

        end if


        ! Fermi level: in case of collinear spin we accept two values (up and down)
        ! call init(fermiBuffer)
        ! call getChildValue(pNode, "FermiLevel", fermiBuffer, errStatus, modifier=modifier)
        ! @:PROPAGATE_ERROR(errStatus)
        ! if ( len(fermiBuffer) .eq. 1) then
        !   call asArray(fermiBuffer, contacts(ii)%eFermi)
        !   contacts(ii)%eFermi(2) = contacts(ii)%eFermi(1)
        ! else if ( len(fermiBuffer) .eq. 2) then
        !   call asArray(fermiBuffer, contacts(ii)%eFermi)
        ! else
        !   call detailedError(pNode, "FermiLevel accepts 1 or 2 (for collinear spin) values",&
        !       & errStatus)
        ! @:PROPAGATE_ERROR(errStatus)
        ! end if
        ! call destruct(fermiBuffer)


        call getChildValue(pNode, "FermiLevel", child1, errStatus, "", child=child2,&
            & allowEmptyValue=.true., modifier=modifier)
        @:PROPAGATE_ERROR(errStatus)
        call getNodeName2(child1, buffer)
        if (char(buffer) == "") then
          contacts(ii)%tFermiSet = .false.
          call detailedWarning(pNode, "Missing Fermi level - required to be set in solver block or&
              & read from a contact shift file")
        else
          call init(fermiBuffer)
          call getChildValue(child2, "", fermiBuffer, errStatus, modifier=modifier)
          @:PROPAGATE_ERROR(errStatus)
          select case(len(fermiBuffer))
          case (1)
            call asArray(fermiBuffer, contacts(ii)%eFermi(1:1))
            contacts(ii)%eFermi(2) = contacts(ii)%eFermi(1)
          case (2)
            call asArray(fermiBuffer, contacts(ii)%eFermi(:2))
          case default
            call detailedError(pNode, "FermiLevel accepts 1 or 2 (for collinear spin) values",&
                & errStatus)
            @:PROPAGATE_ERROR(errStatus)
          end select
          call destruct(fermiBuffer)
          call convertUnitHsd(char(modifier), energyUnits, child2, contacts(ii)%eFermi, errStatus)
          @:PROPAGATE_ERROR(errStatus)

          contacts(ii)%tFermiSet = .true.

          ! NOTE: These options have been commented out: there is a problem in parallel execution
          ! since one single file is accessed by all processors causing rush conditions
          ! The options are therefore disabled for the official dftb+ release
          ! call getChildValue(pNode, "WriteSelfEnergy", contacts(ii)%tWriteSelfEnergy, errStatus,&
          !     & .false.)
          ! @:PROPAGATE_ERROR(errStatus)
          ! call getChildValue(pNode, "WriteSurfaceGF", contacts(ii)%tWriteSurfaceGF, errStatus,&
          !     & .false.)
          ! @:PROPAGATE_ERROR(errStatus)
          ! call getChildValue(pNode, "ReadSelfEnergy", contacts(ii)%tReadSelfEnergy, errStatus,&
          !     & .false.)
          ! @:PROPAGATE_ERROR(errStatus)
          ! call getChildValue(pNode, "ReadSurfaceGF", contacts(ii)%tReadSurfaceGF, errStatus,&
          !     & .false.)
          ! @:PROPAGATE_ERROR(errStatus)
          contacts(ii)%tWriteSelfEnergy = .false.
          contacts(ii)%tWriteSurfaceGF = .false.
          contacts(ii)%tReadSelfEnergy = .false.
          contacts(ii)%tReadSurfaceGF = .false.
        end if

      end if

    end do

  end subroutine readContacts


  !> Read in Fermi levels
  subroutine getFermiLevels(pNode, eFermis, nodeModifier, errStatus)

    !> Document tree node to start from
    type(fnode), pointer :: pNode

    !> Fermi energies for contacts
    real(dp), intent(out) :: eFermis(:)

    !> Any node modifiers in action
    type(string), intent(in) :: nodeModifier

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    real(dp) :: eFermi
    type(fnode), pointer :: pChild
    type(string) :: modifier

    call getChild(pNode, "SetForAll", pChild, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(pChild)) then
      call getChildValue(pChild, "", eFermi, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(nodeModifier), energyUnits, pNode, eFermi, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      eFermis(:) = eFermi
    else
      call getChildValue(pNode, "", eFermis, errStatus, modifier=modifier, child=pChild)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), energyUnits, pChild, eFermis, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine getFermiLevels


  !> Get contacts for terminal currents by name
  subroutine getEmitterCollectorByName(pNode, emitter, collector, contactNames, errStatus)

    !> Node in the input tree for error reporting
    type(fnode), pointer :: pNode

    !> Contact number for emitting
    integer, intent(out) :: emitter

    !> Contact number for collecting
    integer, intent(out) :: collector

    !> Labels of contacts
    character(len=*), intent(in) :: contactNames(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(TListString) :: lString
    character(len=mc) :: buffer

    call init(lString)
    call getChildValue(pNode, "", lString, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    if (len(lString) /= 2) then
      call detailedError(pNode, "You must provide two contacts", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    call get(lString, buffer, 1)
    call getContactByName(contactNames, buffer, pNode, emitter, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call get(lString, buffer, 2)
    call getContactByName(contactNames, buffer, pNode, collector, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call destruct(lString)

  end subroutine getEmitterCollectorByName


  !> Getting the contact by name
  subroutine getContactByName(contactNames, contName, pNode, contact, errStatus)

    !> Node in the input tree for error reporting
    type(fnode), pointer :: pNode

    !> All of the contact labels
    character(len=*), intent(in) :: contactNames(:)

    !> Specific contact label to identify
    character(len=*), intent(in) :: contName

    !> Contact number
    integer, intent(out) :: contact

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    logical :: tFound

    tFound = .false.
    do contact = 1, size(contactNames)
      tFound = (contactNames(contact) == contName)
      if (tFound) then
        exit
      end if
    end do
    if (.not. tFound) then
      call detailedError(pNode, "Invalid collector contact name '" // trim(contName) // "'",&
          & errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine getContactByName


  !> Read the names of regions to calculate PDOS for
  subroutine readPDOSRegions(node, geo, idxdevice, iAtInregion, tShellResInRegion, regionLabels,&
      & errStatus)

    !> Node to be parsed
    type(fnode), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Is the region to be projected by shell
    integer, intent(in) :: idxdevice(2)

    !> Atoms in a given region
    type(TWrappedInt1), allocatable, intent(out) :: iAtInRegion(:)

    !> Is the region to be projected by shell
    logical, allocatable, intent(out) :: tShellResInRegion(:)

    !> Labels for the regions
    character(lc), allocatable, intent(out) :: regionLabels(:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: nReg, iReg
    integer, allocatable :: tmpI1(:)
    type(fnodeList), pointer :: children
    type(fnode), pointer :: child, child2
    type(string) :: buffer
    character(lc) :: strTmp
    logical :: do_ldos

    call getChildren(node, "Region", children)
    nReg = getLength(children)

    if (nReg == 0) then
      call getChildValue(node, "ComputeLDOS", do_ldos, errStatus, .true.)
      @:PROPAGATE_ERROR(errStatus)
      if (do_ldos) then
        write(strTmp,"(I0, ':', I0)") idxdevice(1), idxdevice(2)
        call setChild(node, "Region", child, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call setChildValue(child, "Atoms", trim(strTmp), errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call setChildValue(child, "Label", "localDOS", errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call destroyNodeList(children)
        call getChildren(node, "Region", children)
        nReg = getLength(children)
      else
        return
      end if
    end if

    allocate(tShellResInRegion(nReg))
    allocate(regionLabels(nReg))
    allocate(iAtInRegion(nReg))
    do iReg = 1, nReg
      call getItem1(children, iReg, child)
      call getChildValue(child, "Atoms", buffer, errStatus, child=child2, multiple=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getSelectedAtomIndices(child2, char(buffer), geo%speciesNames,&
          & geo%species(idxdevice(1) : idxdevice(2)), tmpI1, errStatus,&
          & selectionRange=[idxdevice(1), idxdevice(2)], indexRange=[1, geo%nAtom])
      @:PROPAGATE_ERROR(errStatus)
      iAtInRegion(iReg)%data = tmpI1
      call getChildValue(child, "ShellResolved", tShellResInRegion(iReg), errStatus, .false.,&
          & child=child2)
      @:PROPAGATE_ERROR(errStatus)
      if (tShellResInRegion(iReg)) then
        if (.not. all(geo%species(tmpI1) == geo%species(tmpI1(1)))) then
          call detailedError(child2, "Shell resolved PDOS can only summed up over atoms of the same&
              & type", errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end if
      end if
      write(strTmp, "('region',I0)") iReg
      call getChildValue(child, "Label", buffer, errStatus, trim(strTmp))
      @:PROPAGATE_ERROR(errStatus)
      regionLabels(iReg) = unquote(char(buffer))
    end do
    call destroyNodeList(children)

  end subroutine readPDOSRegions


  !> Some assignment and consistency check in negf/poisson containers before calling initialization
  subroutine finalizeNegf(input, errStatus)

    !> Input structure for DFTB+
    type(TInputData), intent(inout) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: ii

    !! Check consistency between different deltas
    if (input%ginfo%tundos%defined.and.input%ginfo%greendens%defined) then
      if (input%ginfo%tundos%delta.ne.input%ginfo%greendens%delta) then
        @:RAISE_ERROR(errStatus, -1, "Delta parameter must be the same in GreensFunction and&
            & TunnelingAndDos")
      end if
    end if

    !! Assign spin degeneracy to every block which may use it
    if (input%ginfo%tundos%defined) then
      if (input%ctrl%tSpin) input%ginfo%tundos%gSpin = 1
      if (.not.input%ctrl%tSpin) input%ginfo%tundos%gSpin = 2
    end if
    if (input%ginfo%greendens%defined) then
      if (input%ctrl%tSpin) input%ginfo%greendens%gSpin = 1
      if (.not.input%ctrl%tSpin) input%ginfo%greendens%gSpin = 2
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Inheritance of first layer indexes to green solver when transport is defined
    if (input%transpar%defined .and. input%ginfo%greendens%defined) then
      input%ginfo%greendens%nPLs = input%transpar%nPLs
      input%ginfo%greendens%PL = input%transpar%PL
    end if

    #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
      !! Not orthogonal directions in transport are only allowed if no Poisson
      if (input%poisson%defined.and.input%transpar%defined) then
        do ii = 1, input%transpar%ncont
          ! If dir is  any value but x,y,z (1,2,3) it is considered oriented along
          ! a direction not parallel to any coordinate axis
          if (input%transpar%contacts(ii)%dir.lt.1 .or. &
            &input%transpar%contacts(ii)%dir.gt.3 ) then
            @:RAISE_ERROR(errStatus, -1, "Contact " // i2c(ii) // " not parallel to any &
              & coordinate axis is not compatible with Poisson solver")
          end if
        end do
      end if
    #:endblock

    !! Temporarily not supporting surface green function read/load
    !! for spin polarised, because spin is handled outside of libnegf
    if (input%ginfo%greendens%defined) then
      if (input%ctrl%tSpin .and. input%ginfo%greendens%saveSGF) then
        @:RAISE_ERROR(errStatus, -1, "SaveSurfaceGFs must be disabled in collinear spin&
            & calculations")
      end if
      if  (input%ctrl%tSpin .and. input%ginfo%greendens%readSGF) then
        @:RAISE_ERROR(errStatus, -1, "ReadSurfaceGFs must be disabled in collinear spin&
            & calculations")
      end if
    end if

  end subroutine finalizeNegf
#:endif


  !> This subroutine overrides the neutral (reference) atom electronic occupation
  subroutine readCustomReferenceOcc(root, orb, referenceOcc, geo, iAtInRegion, customOcc, errStatus)

    !> Node to be parsed
    type(fnode), pointer, intent(in) :: root

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Default reference occupations
    real(dp), intent(in) :: referenceOcc(:,:)

    !> Geometry information
    type(TGeometry), intent(in) :: geo

    !> Atom indices corresponding to user defined reference atomic charges
    type(TWrappedInt1), allocatable, intent(out) :: iAtInRegion(:)

    !> User-defined reference atomic charges
    real(dp), allocatable, intent(out) :: customOcc(:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: node, container, child
    type(fnodeList), pointer :: nodes
    type(string) :: buffer
    integer :: nCustomOcc, iCustomOcc, iShell, iSpecies, nAtom
    character(sc), allocatable :: shellNamesTmp(:)
    logical, allocatable :: atomOverriden(:)

    call renameChildren(root, "CustomizedOccupations", "CustomisedOccupations")
    call getChild(root, "CustomisedOccupations", container, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (.not. associated(container)) then
      return
    end if

    call getChildren(container, "ReferenceOccupation", nodes)
    nCustomOcc = getLength(nodes)
    nAtom = size(geo%species)
    allocate(iAtInRegion(nCustomOcc))
    allocate(customOcc(orb%mShell, nCustomOcc))
    allocate(atomOverriden(nAtom))
    atomOverriden(:) = .false.
    customOcc(:,:) = 0.0_dp

    do iCustomOcc = 1, nCustomOcc
      call getItem1(nodes, iCustomOcc, node)
      call getChildValue(node, "Atoms", buffer, errStatus, child=child, multiple=.true.)
      @:PROPAGATE_ERROR(errStatus)
      call getSelectedAtomIndices(child, char(buffer), geo%speciesNames, geo%species,&
          & iAtInRegion(iCustomOcc)%data, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      if (any(atomOverriden(iAtInRegion(iCustomOcc)%data))) then
        call detailedError(child, "Atom region contains atom(s) which have already been&
            & overridden", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      atomOverriden(iAtInRegion(iCustomOcc)%data) = .true.
      iSpecies = geo%species(iAtInRegion(iCustomOcc)%data(1))
      if (any(geo%species(iAtInRegion(iCustomOcc)%data) /= iSpecies)) then
        call detailedError(child, "All atoms in a ReferenceOccupation declaration must have the&
            & same type.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getShellNames(iSpecies, orb, shellNamesTmp)
      do iShell = 1, orb%nShell(iSpecies)
        call getChildValue(node, shellNamesTmp(iShell), customOcc(iShell, iCustomOcc),&
            & errStatus, default=referenceOcc(iShell, iSpecies))
          @:PROPAGATE_ERROR(errStatus)
      end do
      deallocate(shellNamesTmp)
    end do
    call destroyNodeList(nodes)

  end subroutine readCustomReferenceOcc


  !> Reads the parallel block.
  subroutine readParallel(root, input, errStatus)

    !> Root node eventually containing the current block
    type(fnode), pointer, intent(in) :: root

    !> Input structure to be filled
    type(TInputData), intent(inout) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: node, pTmp

    call getChild(root, "Parallel", node, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (withMpi .and. .not. associated(node)) then
      call setChild(root, "Parallel", node, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    if (associated(node)) then
      if (.not. withMpi) then
        call detailedWarning(node, "Settings will be read but ignored (compiled without MPI&
            & support)")
      end if
      allocate(input%ctrl%parallelOpts)
      call getChildValue(node, "Groups", input%ctrl%parallelOpts%nGroup, errStatus, 1, child=pTmp)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(node, "UseOmpThreads", input%ctrl%parallelOpts%tOmpThreads, errStatus,&
          & .not. withMpi)
      @:PROPAGATE_ERROR(errStatus)
      call readBlacs(node, input%ctrl%parallelOpts%blacsOpts, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine readParallel


  !> Reads the blacs block
  subroutine readBlacs(root, blacsOpts, errStatus)

    !> Root node eventually containing the current block
    type(fnode), pointer, intent(in) :: root

    !> Blacs settings
    type(TBlacsOpts), intent(inout) :: blacsOpts

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: node

    call getChild(root, "Blacs", node, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (withScalapack .and. .not. associated(node)) then
      call setChild(root, "Blacs", node, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    if (associated(node)) then
      if (.not. withScalapack) then
        call detailedWarning(node, "Settings will be read but ignored (compiled without SCALAPACK&
            & support)")
      end if
      call getChildValue(node, "BlockSize", blacsOpts%blockSize, errStatus, 32)
      @:PROPAGATE_ERROR(errStatus)
    end if

  end subroutine readBlacs


  !> Reads the settings for electrostatic potential plotting
  subroutine readElectrostaticPotential(node, geo, ctrl, errStatus)

    !> Node containing optional electrostatic settings
    type(fnode), pointer, intent(in) :: node

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Control structure
    type(TControl), intent(inout) :: ctrl

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child, child2, child3
    type(string) :: buffer, modifier
    type(TListRealR1) :: lr1

    call getChild(node, "ElectrostaticPotential", child, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (.not. associated(child)) then
      return
    end if

    if (.not. ctrl%tSCC) then
      @:RAISE_ERROR(errStatus, -1, "Electrostatic potentials only available in an SCC calculation")
    end if
    allocate(ctrl%elStatPotentialsInp)
    call getChildValue(child, "OutputFile", buffer, errStatus, "ESP.dat")
    @:PROPAGATE_ERROR(errStatus)
    ctrl%elStatPotentialsInp%espOutFile = unquote(char(buffer))
    ctrl%elStatPotentialsInp%tAppendEsp = .false.
    if (ctrl%isGeoOpt .or. ctrl%tMD) then
      call getChildValue(child, "AppendFile", ctrl%elStatPotentialsInp%tAppendEsp, errStatus,&
          & .false.)
      @:PROPAGATE_ERROR(errStatus)
    end if
    call init(lr1)
    ! discrete points
    call getChildValue(child, "Points", child2, errStatus, "", child=child3, modifier=modifier,&
        & allowEmptyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName2(child2, buffer)
    if (char(buffer) /= "") then
      call getChildValue(child3, "", 3, lr1, errStatus, modifier=modifier)
      @:PROPAGATE_ERROR(errStatus)
      allocate(ctrl%elStatPotentialsInp%espGrid(3,len(lr1)))
      call asArray(lr1, ctrl%elStatPotentialsInp%espGrid)
      if (geo%tPeriodic .and. (char(modifier) == "F" .or. char(modifier) == "f")) then
        ctrl%elStatPotentialsInp%espGrid = matmul(geo%latVecs, ctrl%elStatPotentialsInp%espGrid)
      else
        call convertUnitHsd(char(modifier), lengthUnits, child3, ctrl%elStatPotentialsInp%espGrid,&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end if
    call destruct(lr1)

    ! grid specification for points instead
    call getChild(child, "Grid", child2, errStatus, modifier=modifier, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (associated(child2)) then
      if (allocated(ctrl%elStatPotentialsInp%espGrid)) then
        @:RAISE_ERROR(errStatus, -1, "Both grid and point specification not both currently&
            & possible")
      end if
      if (geo%tPeriodic) then
        call readGrid(ctrl%elStatPotentialsInp%espGrid, child2, modifier, errStatus,&
            & latVecs=geo%latVecs, nPoints=ctrl%elStatPotentialsInp%gridDimensioning,&
            & origin=ctrl%elStatPotentialsInp%origin,&
            & axes=ctrl%elStatPotentialsInp%axes)
        @:PROPAGATE_ERROR(errStatus)
      else
        call readGrid(ctrl%elStatPotentialsInp%espGrid, child2, modifier, errStatus,&
            & nPoints=ctrl%elStatPotentialsInp%gridDimensioning,&
            & origin=ctrl%elStatPotentialsInp%origin,&
            & axes=ctrl%elStatPotentialsInp%axes)
        @:PROPAGATE_ERROR(errStatus)
      end if
    end if
    if (.not.allocated(ctrl%elStatPotentialsInp%espGrid)) then
      call detailedError(child,"Either a grid or set of points must be specified", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    call getChildValue(child, "Softening", ctrl%elStatPotentialsInp%softenESP, errStatus,&
        & 1.0E-6_dp, modifier=modifier, child=child2)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), lengthUnits, child2, ctrl%elStatPotentialsInp%softenEsp,&
        & errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readElectrostaticPotential


  !> Read in a grid specification
  subroutine readGrid(points, node, modifier, errStatus, latVecs, nPoints, origin, axes)

    !> Points in the grid
    real(dp), allocatable, intent(out) :: points(:,:)

    !> input data to parse
    type(fnode), pointer, intent(in) :: node

    !> unit modifier for the grid
    type(string), intent(in) :: modifier

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !> geometry of the system
    real(dp), intent(in), optional :: latVecs(:,:)

    !> Number of grid points in each direction, if required
    integer, intent(out), optional :: nPoints(3)

    !> origin of grid if required
    real(dp), intent(out), optional :: origin(3)

    !> axes of the grid if required
    real(dp), intent(out), optional :: axes(3,3)

    type(fnode), pointer :: child
    real(dp) :: r3Tmp(3), r3Tmpb(3)
    integer :: i3Tmp(3), iPt, ii, jj, kk
    logical :: tPeriodic
    real(dp) :: axes_(3,3), r33Tmp(3,3)

    tPeriodic = present(latvecs)

    if (.not.tPeriodic .and. (char(modifier) == "F" .or. char(modifier) == "f")) then
      call detailedError(node, "Fractional grid specification only available for periodic&
          & geometries", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    call getChildValue(node, "Spacing", r3Tmp, errStatus, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "Origin", r3Tmpb, errStatus, child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "GridPoints", i3Tmp, errStatus, child=child)
    @:PROPAGATE_ERROR(errStatus)
    if (any(i3Tmp < 1)) then
      call detailedError(child, "Grid must be at least 1x1x1", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    if (any(abs(r3Tmp) < epsilon(1.0_dp) .and. i3Tmp > 1)) then
      call detailedError(child, "Grid spacings must be non-zero", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if
    allocate(points(3,product(i3Tmp)))
    if (present(nPoints)) then
      nPoints = i3Tmp
    end if

    !  length not fraction modifier
    if (.not.(tPeriodic .and. (char(modifier) == "F" .or. char(modifier) == "f"))) then
      call convertUnitHsd(char(modifier), lengthUnits, child, r3Tmp, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call convertUnitHsd(char(modifier), lengthUnits, child, r3Tmpb, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    points = 0.0_dp
    iPt = 0
    do ii = 0, i3Tmp(1)-1
      do jj = 0, i3Tmp(2)-1
        do kk = 0, i3Tmp(3)-1
          iPt = iPt + 1
          points(1,iPt) = ii * r3Tmp(1) + r3Tmpb(1)
          points(2,iPt) = jj * r3Tmp(2) + r3Tmpb(2)
          points(3,iPt) = kk * r3Tmp(3) + r3Tmpb(3)
        end do
      end do
    end do

    ! transformation matrix on directions, could use a 4x4 homogeneous coordinate transform instead
    if (.not.(char(modifier) == "F" .or. char(modifier) == "f") .or. .not.tPeriodic) then
      r33Tmp = reshape([1,0,0,0,1,0,0,0,1],[3,3])
      call getChildValue(node, "Directions", axes_, errStatus, r33Tmp, child=child)
      @:PROPAGATE_ERROR(errStatus)
      if (abs(determinant33(axes_)) < epsilon(1.0_dp)) then
        call detailedError(child, "Dependent axis directions", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      do ii = 1, 3
        axes_(:,ii) = axes_(:,ii) / sqrt(sum(axes_(:,ii)**2))
      end do
      points = matmul(axes_,points)
      if (present(axes)) then
        axes = axes_*spread(r3Tmp,2,3)
      end if
    end if

    if (present(origin)) then
      origin = r3Tmpb
    end if

    ! Fractional specification of points
    if (tPeriodic .and. (char(modifier) == "F" .or. char(modifier) == "f")) then
      points = matmul(latVecs,points)
      if (present(origin)) then
        origin = matmul(latVecs,origin)
      end if
      if (present(axes)) then
        axes = latVecs * spread(r3Tmp,2,3)
      end if
    end if

  end subroutine readGrid

  function is_numeric(string) result(is)
    character(len=*), intent(in) :: string
    logical :: is

    real :: x
    integer :: err

    read(string,*,iostat=err) x
    is = (err == 0)
  end function is_numeric


  !> Parse range separation input
  subroutine parseRangeSeparated(node, input, skFiles, errStatus)
    !> Node to parse
    type(fnode), intent(in), pointer :: node

    !> Range separated data structure to fill
    type(TRangeSepInp), intent(inout), allocatable :: input

    !> List of SK file names to read in for every interaction
    type(TListCharLc), intent(inout) :: skFiles(:,:)

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    !! File name of representative SK-file to read
    character(lc) :: fileName

    !! Range-separated extra tag in SK-files, if allocated
    integer :: rangeSepSkTag

    !! Range-separated functional type of user input
    integer :: rangeSepInputTag

    !! Auxiliary node pointers
    type(fnode), pointer :: hybridChild, hybridValue, screeningChild, screeningValue, child1

    !! True, if RangeSeparated input block present
    logical :: isHybridInp

    !! True, if range-separated extra tag found in SK-file(s)
    logical :: isHybridSk

    !! Temporary string buffers
    type(string) :: buffer, modifier

    @:ASSERT(size(skFiles, dim=1) == size(skFiles, dim=2))
    @:ASSERT((size(skFiles, dim=1) > 0))

    ! Extracting range-separated tag from first SK-file only is a workaround and assumes that a set
    ! of given SK-files uses the same parameters (which should always be the case)!
    call get(skFiles(1, 1), fileName, 1)

    ! Check if SK-files contain extra tag for range-separated functionals
    call inquireRangeSepTag(fileName, rangeSepSkTag)
    isHybridSk = rangeSepSkTag /= rangeSepFunc%none

    call getChild(node, "RangeSeparated", hybridChild, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    isHybridInp = associated(hybridChild)

    if (isHybridInp .and. .not. isHybridSk) then
      @:RAISE_ERROR(errStatus, -1, "RangeSeparated input block present, but SK-file '"&
          & // trim(fileName) // "' seems to be (semi-)local.")
    elseif (isHybridSk .and. .not. isHybridInp) then
      @:RAISE_ERROR(errStatus, -1, "Hybrid SK-file '" // trim(fileName) // "' present, but HSD&
          & input block missing.")
    end if

    if (isHybridInp) then
      call getChildValue(node, "RangeSeparated", hybridValue, errStatus, "None", child=hybridChild)
      @:PROPAGATE_ERROR(errStatus)
      call getNodeName(hybridValue, buffer)
      ! Convert hybrid functional type of user input to enumerator
      select case(tolower(char(buffer)))
      case ("lc")
        rangeSepInputTag = rangeSepFunc%lc
      case default
        call detailedError(hybridChild, "Unknown hybrid xc-functional type '" // char(buffer)&
            & // "' in input.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end select
      if (.not. rangeSepInputTag == rangeSepSkTag) then
        ! Check if hybrid functional type is in line with SK-files
        call detailedError(hybridChild, "Hybrid functional type conflict with SK-files.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      select case (tolower(char(buffer)))
      case ("none")
        rangeSepInputTag = rangeSepFunc%none
        continue
      case ("lc")
        allocate(input)
        call getChildValue(hybridValue, "Screening", screeningValue, errStatus, "Thresholded",&
            & child=screeningChild)
        @:PROPAGATE_ERROR(errStatus)
        call getNodeName(screeningValue, buffer)
        select case(char(buffer))
        case ("neighbourbased")
          input%rangeSepAlg = rangeSepTypes%neighbour
          call getChildValue(screeningValue, "CutoffReduction", input%cutoffRed, errStatus, 0.0_dp,&
              & modifier=modifier, child=child1)
          @:PROPAGATE_ERROR(errStatus)
          call convertUnitHsd(char(modifier), lengthUnits, child1, input%cutoffRed, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        case ("thresholded")
          input%rangeSepAlg = rangeSepTypes%threshold
          call getChildValue(screeningValue, "Threshold", input%screeningThreshold, errStatus,&
              & 1e-6_dp)
          @:PROPAGATE_ERROR(errStatus)
          call getChildValue(screeningValue, "CutoffReduction", input%cutoffRed, errStatus, 0.0_dp,&
              & modifier=modifier, child=child1)
          @:PROPAGATE_ERROR(errStatus)
          call convertUnitHsd(char(modifier), lengthUnits, child1, input%cutoffRed, errStatus)
          @:PROPAGATE_ERROR(errStatus)
        case ("matrixbased")
          input%rangeSepAlg = rangeSepTypes%matrixBased
          ! In this case, CutoffRedunction is not used so it should be set to zero.
          input%cutoffRed = 0.0_dp
        case default
          call getNodeHSdName(screeningValue, buffer)
          call detailedError(screeningChild, "Invalid screening method '" // char(buffer) // "'.",&
              & errStatus)
          @:PROPAGATE_ERROR(errStatus)
        end select
      end select
      input%rangeSepType = rangeSepInputTag
    end if

  end subroutine parseRangeSeparated


  !> Reads the REKS block
  subroutine readReks(node, dummy, ctrl, geo, errStatus)

    !> Node to parse
    type(fnode), pointer, intent(in) :: node

    !> Node to parse
    type(fnode), pointer, intent(in) :: dummy

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(string) :: buffer

    ! SSR(2,2) or SSR(4,4) stuff
    call getNodeName(dummy, buffer)

    select case (char(buffer))
    case ("none")
      ctrl%reksInp%reksAlg = reksTypes%noReks
    case ("ssr22")
      ctrl%reksInp%reksAlg = reksTypes%ssr22
      call readSSR22(dummy, ctrl, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case ("ssr44")
      ctrl%reksInp%reksAlg = reksTypes%ssr44
      call detailedError(node, "SSR(4,4) is not implemented yet.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case default
      call getNodeHSDName(dummy, buffer)
      call detailedError(node, "Invalid Algorithm '" // char(buffer) // "'", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

  end subroutine readReks


  !> Reads the SSR(2,2) block
  subroutine readSSR22(node, ctrl, geo, errStatus)

    !> Node to parse
    type(fnode), pointer, intent(in) :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child1, value2, child2
    type(TListString) :: strBuffer
    type(string) :: buffer2
    character(sc), allocatable :: tmpFunc(:)
    integer :: ii, nFunc
    logical :: tFunc = .true.


    !> Read 'Energy' block
    call getChild(node, "Energy", child1, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    !> Read 'Functional' block in 'Energy' block
    call init(strBuffer)
    call getChildValue(child1, "Functional", strBuffer, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    allocate(tmpFunc(len(strBuffer)))
    call asArray(strBuffer, tmpFunc)
    call destruct(strBuffer)

    !> Decide the energy functionals to be included in SA-REKS(2,2)
    nFunc = size(tmpFunc, dim=1)
    if (nFunc == 1) then
      if (trim(tmpFunc(1)) == "PPS") then
        !> Minimized energy functional : PPS
        ctrl%reksInp%Efunction = 1
      else
        tFunc = .false.
      end if
    else if (nFunc == 2) then
      if (trim(tmpFunc(1)) == "PPS" .and. trim(tmpFunc(2)) == "OSS") then
        !> Minimized energy functional : (PPS+OSS)/2
        ctrl%reksInp%Efunction = 2
      else
        tFunc = .false.
      end if
    else
      tFunc = .false.
    end if

    if (.not. tFunc) then
      write(stdOut,'(A)',advance="no") "Current Functional : "
      do ii = 1, nFunc
        if (ii == nFunc) then
          write(stdOut,'(A)') "'" // trim(tmpFunc(ii)) // "'"
        else
          write(stdOut,'(A)',advance="no") "'" // trim(tmpFunc(ii)) // "' "
        end if
      end do
      call detailedError(child1, "Invalid Functional", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end if

    !> Decide the energy states in SA-REKS
    !> If true, it includes all possible states in current active space
    !> If false, it includes the states used in minimized energy functional
    call getChildValue(child1, "IncludeAllStates", ctrl%reksInp%tAllStates, errStatus,&
        & default=.false.)
    @:PROPAGATE_ERROR(errStatus)
    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    call getChildValue(child1, "StateInteractions", ctrl%reksInp%tSSR, errStatus, default=.false.)
    @:PROPAGATE_ERROR(errStatus)

    !> Target SSR state
    call getChildValue(node, "TargetState", ctrl%reksInp%rstate, errStatus, default=1)
    @:PROPAGATE_ERROR(errStatus)
    !> Target microstate
    call getChildValue(node, "TargetMicrostate", ctrl%reksInp%Lstate, errStatus, default=0)
    @:PROPAGATE_ERROR(errStatus)

    !> Read initial guess for eigenvectors in REKS
    !> If true, initial eigenvectors are obtained from 'eigenvec.bin'
    !> If false, initial eigenvectors are obtained from diagonalisation of H0
    call getChildValue(node, "ReadEigenvectors", ctrl%reksInp%tReadMO, errStatus, default=.false.)
    @:PROPAGATE_ERROR(errStatus)
    !> Maximum iteration used in FON optimisation
    call getChildValue(node, "FonMaxIter", ctrl%reksInp%FonMaxIter, errStatus, default=20)
    @:PROPAGATE_ERROR(errStatus)
    !> Shift value in SCC cycle
    call getChildValue(node, "Shift", ctrl%reksInp%shift, errStatus, default=0.3_dp)
    @:PROPAGATE_ERROR(errStatus)

    !> Read "SpinTuning" block with 'nType' elements
    call readSpinTuning(node, ctrl, geo%nSpecies, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    !> Calculate transition dipole moments
    call getChildValue(node, "TransitionDipole", ctrl%reksInp%tTDP, errStatus, default=.false.)
    @:PROPAGATE_ERROR(errStatus)


    !> Read 'Gradient' block
    !> Algorithms to calculate analytical gradients
    call getChildValue(node, "Gradient", value2, errStatus, "ConjugateGradient", child=child2)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName(value2, buffer2)

    select case (char(buffer2))
    case ("conjugategradient")
      !> Maximum iteration used in calculation of gradient with PCG and CG
      call getChildValue(value2, "CGmaxIter", ctrl%reksInp%CGmaxIter, errStatus, default=20)
      @:PROPAGATE_ERROR(errStatus)
      !> Tolerance used in calculation of gradient with PCG and CG
      call getChildValue(value2, "Tolerance", ctrl%reksInp%Glimit, errStatus, default=1.0E-8_dp)
      @:PROPAGATE_ERROR(errStatus)
      !> Use preconditioner for conjugate gradient algorithm
      call getChildValue(value2, "Preconditioner", ctrl%reksInp%tPrecond, errStatus,&
          & default=.false.)
      @:PROPAGATE_ERROR(errStatus)
      !> Save 'A' and 'Hxc' to memory in gradient calculation
      call getChildValue(value2, "SaveMemory", ctrl%reksInp%tSaveMem, errStatus, default=.false.)
      @:PROPAGATE_ERROR(errStatus)
      if (ctrl%reksInp%tPrecond) then
        !> 1: preconditioned conjugate gradient (PCG)
        ctrl%reksInp%Glevel = 1
      else
        !> 2: conjugate gradient (CG)
        ctrl%reksInp%Glevel = 2
      end if
    case ("direct")
      !> 3: direct inverse-matrix multiplication
      ctrl%reksInp%Glevel = 3
    case default
      call getNodeHSDName(value2, buffer2)
      call detailedError(child2, "Invalid Algorithm '" // char(buffer2) // "'", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    end select

    !> Calculate relaxed density of SSR or SA-REKS state
    call getChildValue(node, "RelaxedDensity", ctrl%reksInp%tRD, errStatus, default=.false.)
    @:PROPAGATE_ERROR(errStatus)
    !> Calculate nonadiabatic coupling vectors
    call getChildValue(node, "NonAdiabaticCoupling", ctrl%reksInp%tNAC, errStatus, default=.false.)
    @:PROPAGATE_ERROR(errStatus)

    !> Print level in standard output file
    call getChildValue(node, "VerbosityLevel", ctrl%reksInp%Plevel, errStatus, default=1)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readSSR22


  !> Reads SpinTuning block in REKS input
  subroutine readSpinTuning(node, ctrl, nType, errStatus)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Number of types for atoms
    integer, intent(in) :: nType

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: value1, child
    type(string) :: buffer, modifier
    type(TListRealR1) :: realBuffer
    integer :: nAtom, iType
    real(dp), allocatable :: tmpTuning(:,:)

    call getChildValue(node, "SpinTuning", value1, errStatus, "", child=child, modifier=modifier,&
        & allowEmptyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call getNodeName2(value1, buffer)
    if (char(buffer) == "") then
      ! no 'SpinTuning' block in REKS input
      allocate(ctrl%reksInp%Tuning(nType))
      do iType = 1, nType
        ctrl%reksInp%Tuning(iType) = 1.0_dp
      end do
    else
      ! 'SpinTuning' block in REKS input
      call init(realBuffer)
      call getChildValue(child, "", 1, realBuffer, errStatus, modifier=modifier)
      @:PROPAGATE_ERROR(errStatus)
      nAtom = len(realBuffer)
      if (nAtom /= nType) then
        call detailedError(node, "Incorrect number of 'SpinTuning' block: "&
            & // i2c(nAtom) // " supplied, "&
            & // i2c(nType) // " required.", errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      allocate(tmpTuning(1,nAtom))
      call asArray(realBuffer, tmpTuning)
      call destruct(realBuffer)
      allocate(ctrl%reksInp%Tuning(nType))
      ctrl%reksInp%Tuning(:) = tmpTuning(1,:)
    end if

  end subroutine readSpinTuning


  !> Parses Chimes related options.
  subroutine parseChimes(root, chimesRepInput, errStatus)
    type(fnode), pointer, intent(in) :: root
    type(TChimesRepInp), allocatable, intent(out) :: chimesRepInput

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: chimes
    type(string) :: buffer

  #:if WITH_CHIMES
    type(string), allocatable :: searchPath(:)
  #:endif

    character(len=:), allocatable :: chimesFile

    call getChild(root, "Chimes", chimes, errStatus, requested=.false.)
    @:PROPAGATE_ERROR(errStatus)
    if (.not. associated(chimes)) return
    #:if WITH_CHIMES
      allocate(chimesRepInput)
      call getChildValue(chimes, "ParameterFile", buffer, errStatus, default="chimes.dat")
      @:PROPAGATE_ERROR(errStatus)
      chimesFile = unquote(char(buffer))
      call getParamSearchPath(searchPath)
      call findFile(searchPath, chimesFile, chimesRepInput%chimesFile)
      if (.not. allocated(chimesRepInput%chimesFile)) then
        @:RAISE_ERROR(errStatus, -1, "Could not find ChIMES parameter file '" // chimesFile // "'")
      end if
    #:else
      call detailedError(chimes, "ChIMES repuslive correction requested, but code was compiled&
          & without ChIMES support", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    #:endif

  end subroutine parseChimes


  !> Returns parser version for a given input version or throws an error if not possible.
  subroutine parserVersionFromInputVersion(versionString, node, parserVersion, errStatus)

    !> Input version string
    character(len=*), intent(in) :: versionString

    !> Input version node (needed for error messagess)
    type(fnode), pointer :: node

    !> Corresponding parser version.
    integer, intent(out) :: parserVersion

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    integer :: ii

    do ii = 1, size(versionMaps)
      if (versionMaps(ii)%inputVersion == versionString) then
        parserVersion = versionMaps(ii)%parserVersion
        return
      end if
    end do

    call detailedError(node, "Program version '"// trim(versionString) // "' is not recognized",&
        & errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine parserVersionFromInputVersion


end module dftbp_dftbplus_parser
