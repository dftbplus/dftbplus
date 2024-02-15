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
  use dftbp_dftb_periodic, only : TNeighbourList, TNeighbourlist_init, getSuperSampling,&
      & getCellTranslations, updateNeighbourList
  use dftbp_dftb_hybridxc, only : THybridXcSKTag, hybridXcAlgo, hybridXcGammaTypes,&
      & checkSupercellFoldingMatrix, hybridXcFunc
  use dftbp_dftb_repulsive_chimesrep, only : TChimesRepInp
  use dftbp_dftb_repulsive_polyrep, only : TPolyRepInp, TPolyRep
  use dftbp_dftb_repulsive_splinerep, only : TSplineRepInp, TSplineRep
  use dftbp_dftb_slakocont, only : init, addTable
  use dftbp_dftb_slakoeqgrid, only : skEqGridNew, skEqGridOld, TSlakoEqGrid, init
  use dftbp_dftbplus_forcetypes, only : forceTypes
  use dftbp_dftbplus_input_fileaccess, only : readBinaryAccessTypes
  use dftbp_dftbplus_inputconversion, only : transformpdosregioninfo
  use dftbp_dftbplus_inputdata, only :TInputData, TControl, TSlater, TBlacsOpts, THybridXcInp
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
  use dftbp_type_oldskdata, only : TOldSKData, readFromFile, inquireHybridXcTag
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
  subroutine readHsdFile(hsdFile, hsdTree)

    !> Name of the input file
    character(*), intent(in) :: hsdFile

    !> Data tree representation of the input
    type(fnode), pointer :: hsdTree

    call parseHSD(rootTag, hsdFile, hsdTree)

  end subroutine readHsdFile


  !> Parse input from an HSD/XML file
  subroutine parseHsdTree(hsdTree, input, parserFlags)

    !> Tree representation of the input
    type(fnode), pointer :: hsdTree

    !> Returns initialised input variables on exit
    type(TInputData), intent(out) :: input

    !> Special block containings parser related settings
    type(TParserFlags), intent(out) :: parserFlags

    type(TStatus) :: errStatus
    type(TOrbitals) :: orb
    type(fnode), pointer :: root, tmp, driverNode, hamNode, analysisNode, child, dummy
    logical :: tReadAnalysis
    integer, allocatable :: implicitParserVersion

    write(stdout, '(A,1X,I0,/)') 'Parser version:', parserVersion
    write(stdout, "(A)") repeat("-", 80)

    call getChild(hsdTree, rootTag, root)

    call handleInputVersion(root, implicitParserVersion)
    call getChildValue(root, "ParserOptions", dummy, "", child=child, list=.true.,&
        & allowEmptyValue=.true., dummyValue=.true.)
    call readParserOptions(child, root, parserFlags, implicitParserVersion)

    ! Read the geometry unless the list of atoms has been provided through the API
    if (.not. allocated(input%geom%coords)) then
      call getChild(root, "Geometry", tmp)
      call readGeometry(tmp, input)
    end if

    ! Hamiltonian settings that need to know settings from the REKS block
    call getChildValue(root, "Reks", dummy, "None", child=child)
    call readReks(child, dummy, input%ctrl, input%geom)

    call getChild(root, "Transport", child, requested=.false.)

  #:if WITH_TRANSPORT

    ! Read in transport and modify geometry if it is only a contact calculation
    if (associated(child)) then
      call readTransportGeometry(child, input%geom, input%transpar)
    else
      input%transpar%ncont=0
      allocate(input%transpar%contacts(0))
      ! set range of atoms in the 'device', as there are no contacts
      input%transpar%idxdevice(1) = 1
      input%transpar%idxdevice(2) = input%geom%nAtom
    end if

    call getChild(root, "Dephasing", child, requested=.false.)
    if (associated(child)) then
      call detailedError(child, "Be patient... Dephasing feature will be available soon!")
      !call readDephasing(child, input%slako%orb, input%geom, input%transpar, input%ginfo%tundos)
    end if

    ! electronic Hamiltonian
    call getChildValue(root, "Hamiltonian", hamNode)
    call readHamiltonian(hamNode, input%ctrl, input%geom, input%slako, input%transpar,&
        & input%ginfo%greendens, input%poisson, errStatus)

  #:else

    if (associated(child)) then
      call detailedError(child, "Program has been compiled without transport enabled")
    end if

    ! electronic Hamiltonian
    call getChildValue(root, "Hamiltonian", hamNode)
    call readHamiltonian(hamNode, input%ctrl, input%geom, input%slako, input%poisson, errStatus)

  #:endif

    if (errStatus%hasError()) then
      call error(errStatus%message)
    end if

    call getChildValue(root, "Driver", driverNode, "", child=child, allowEmptyValue=.true.)
  #:if WITH_TRANSPORT
    call readDriver(driverNode, child, input%geom, input%ctrl, input%transpar)
  #:else
    call readDriver(driverNode, child, input%geom, input%ctrl)
  #:endif

    tReadAnalysis = .true.
    call getChild(root, "ElectronDynamics", child=child, requested=.false.)
    if (associated(child)) then
      allocate(input%ctrl%elecDynInp)
      call readElecDynamics(child, input%ctrl%elecDynInp, input%geom, input%ctrl%masses)
      if (input%ctrl%elecDynInp%tReadRestart .and. .not.input%ctrl%elecDynInp%tPopulations) then
        tReadAnalysis = .false.
      end if

    end if

    if (tReadAnalysis) then
      ! Analysis of properties
      call getChildValue(root, "Analysis", dummy, "", child=analysisNode, list=.true., &
          & allowEmptyValue=.true., dummyValue=.true.)

    #:if WITH_TRANSPORT
      call readAnalysis(analysisNode, input%ctrl, input%geom, input%slako%orb, input%transpar, &
          & input%ginfo%tundos)

      call finalizeNegf(input)
    #:else
      call readAnalysis(analysisNode, input%ctrl, input%geom)
    #:endif

    end if

    call getChildValue(root, "ExcitedState", dummy, "", child=child, list=.true., &
        & allowEmptyValue=.true., dummyValue=.true.)
    call readExcited(child, input%geom, input%ctrl)

    ! Hamiltonian settings that need to know about settings from the blocks above
    call readLaterHamiltonian(hamNode, input%ctrl, driverNode, input%geom)

    call getChildValue(root, "Options", dummy, "", child=child, list=.true., &
        & allowEmptyValue=.true., dummyValue=.true.)
    call readOptions(child, input%ctrl)

    ! W values if needed by Hamiltonian or excited state calculation
    if (allocated(input%ctrl%tbliteInp)) then
      call input%ctrl%tbliteInp%setupOrbitals(input%geom%species, orb)
      call readSpinConstants(hamNode, input%geom, orb, input%ctrl)
    else
      call readSpinConstants(hamNode, input%geom, input%slako%orb, input%ctrl)
    end if

    ! analysis settings that need to know settings from the options block
    if (tReadAnalysis) then
      call readLaterAnalysis(analysisNode, input%ctrl)
    end if

    ! read parallel calculation settings
    call readParallel(root, input)

    ! input data strucutre has been initialised
    input%tInitialized = .true.

  end subroutine parseHsdTree


  !> Converts input version to parser version and removes InputVersion node if present.
  subroutine handleInputVersion(root, implicitParserVersion)

    !> Root eventually containing InputVersion
    type(fnode), pointer, intent(in) :: root

    !> Parser version corresponding to input version, or unallocated if none has been found
    integer, allocatable, intent(out) :: implicitParserVersion

    type(fnode), pointer :: child, dummy
    type(string) :: versionString

    call getChild(root, "InputVersion", child, requested=.false.)
    if (associated(child)) then
      call getChildValue(child, "", versionString)
      implicitParserVersion = parserVersionFromInputVersion(unquote(char(versionString)), child)
      dummy => removeChild(root, child)
      call destroyNode(dummy)
    end if

  end subroutine handleInputVersion


  !> Read in parser options (options not passed to the main code)
  subroutine readParserOptions(node, root, flags, implicitVersion)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Root of the entire tree (in case it needs to be converted, for example because of
    !> compatibility options)
    type(fnode), pointer :: root

    !> Contains parser flags on exit.
    type(TParserFlags), intent(out) :: flags

    !> Parser version implied by version number
    integer, intent(in), optional :: implicitVersion

    integer :: inputVersion
    type(fnode), pointer :: child

    if (present(implicitVersion)) then
      call getChild(node, "ParserVersion", child, requested=.false.)
      if (associated(child)) then
        call getChildValue(child, "", inputVersion)
        if (inputVersion /= implicitVersion) then
          call detailedError(child, "Parser version deduced from InputVersion ("&
          & // i2c(implicitVersion) // ") differs from version explicitely set in ParserVersion ("&
          & // i2c(inputVersion) // ")")
        end if
      else
        inputVersion = implicitVersion
        call setChildValue(node, "ParserVersion", inputVersion, child=child)
      end if
    else
      call getChildValue(node, "ParserVersion", inputVersion, parserVersion, child=child)
    end if

    if (inputVersion < 1 .or. inputVersion > parserVersion) then
      call detailedError(child, "Invalid parser version (" // i2c(inputVersion) // ")")
    else if (inputVersion < minVersion) then
      call detailedError(child, &
          & "Sorry, no compatibility mode for parser version " // i2c(inputVersion)&
          & // " (too old)")
    else if (inputVersion /= parserVersion) then
      write(stdout, "(A,I2,A,I2,A)") "***  Converting input from parser version ",&
          & inputVersion, " to parser version ", parserVersion, " ..."
      call convertOldHSD(root, inputVersion, parserVersion)
      write(stdout, "(A,/)") "***  Done."
    end if

    call getChildValue(node, "WriteHSDInput", flags%tWriteHSD, .true.)
    if (.not. flags%tWriteHSD) then
      call detailedWarning(node, "WriteHSDInput turned off. You are not guaranteed" // newline // &
          &" to able to obtain the same results with a later version of the code!" // newline // &
          & "(the dftb_pin.hsd file DOES guarantee this)")
    end if
    call getChildValue(node, "StopAfterParsing", flags%tStop, .false.)

    call getChildValue(node, "IgnoreUnprocessedNodes", &
        &flags%tIgnoreUnprocessed, .false.)

  end subroutine readParserOptions


  !> Read in Geometry
  subroutine readGeometry(node, input)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Input structure to be filled
    type(TInputData), intent(inout) :: input

    type(fnode), pointer :: value1, child
    type(string) :: buffer

    call getChildValue(node, "", value1, child=child)
    call getNodeName(value1, buffer)
    input%geom%tPeriodic = .false.
    input%geom%tHelical = .false.
    select case (char(buffer))
    case ("genformat")
      call readTGeometryGen(value1, input%geom)
    case ("xyzformat")
      call readTGeometryXyz(value1, input%geom)
    case ("vaspformat")
      call readTGeometryVasp(value1, input%geom)
    case ("lammpsformat")
      call readTGeometryLammps(value1, input%geom)
    case default
      call setUnprocessed(value1)
      call readTGeometryHSD(child, input%geom)
    end select

  end subroutine readGeometry


  !> Read in driver properties
#:if WITH_TRANSPORT
  subroutine readDriver(node, parent, geom, ctrl, transpar)
#:else
  subroutine readDriver(node, parent, geom, ctrl)
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
            & geometries")
      end if

      allocate(ctrl%geoOpt)

      call readGeoOptInput(node, geom, ctrl%geoOpt, atomsRange)

      call getChildValue(node, "AppendGeometries", ctrl%tAppendGeo, .false.)

      ctrl%tForces = .true.
      ctrl%restartFreq = 1

    case ("steepestdescent")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ! Steepest downhill optimisation
      ctrl%iGeoOpt = geoOptTypes%steepestDesc
      call commonGeoOptions(node, ctrl, geom, atomsRange)

    case ("conjugategradient")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')// "Please use the GeometryOptimisation driver instead.")

      ! Conjugate gradient location optimisation
      ctrl%iGeoOpt = geoOptTypes%conjugateGrad
      call commonGeoOptions(node, ctrl, geom, atomsRange)

    case("gdiis")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ! Gradient DIIS optimisation, only stable in the quadratic region
      ctrl%iGeoOpt = geoOptTypes%diis
      call getChildValue(node, "alpha", ctrl%deltaGeoOpt, 1.0E-1_dp)
      call getChildValue(node, "Generations", ctrl%iGenGeoOpt, 8)
      call commonGeoOptions(node, ctrl, geom, atomsRange)

    case ("lbfgs")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ctrl%iGeoOpt = geoOptTypes%lbfgs

      allocate(ctrl%lbfgsInp)
      call getChildValue(node, "Memory", ctrl%lbfgsInp%memory, 20)

      call getChildValue(node, "LineSearch", ctrl%lbfgsInp%isLineSearch, .false.)

      isMaxStepNeeded = .not. ctrl%lbfgsInp%isLineSearch
      if (isMaxStepNeeded) then
        call getChildValue(node, "setMaxStep", ctrl%lbfgsInp%isLineSearch, isMaxStepNeeded)
        ctrl%lbfgsInp%MaxQNStep = isMaxStepNeeded
      else
        call getChildValue(node, "oldLineSearch", ctrl%lbfgsInp%isOldLS, .false.)
      end if

      call commonGeoOptions(node, ctrl, geom, atomsRange, ctrl%lbfgsInp%isLineSearch)

    case ("fire")

      modeName = "geometry relaxation"
      call detailedWarning(node, "This driver is deprecated and will be removed in future&
          & versions."//new_line('a')//&
          & "Please use the GeometryOptimisation driver instead.")

      ctrl%iGeoOpt = geoOptTypes%fire
      call commonGeoOptions(node, ctrl, geom, atomsRange, .false.)
      call getChildValue(node, "TimeStep", ctrl%deltaT, 1.0_dp, modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), timeUnits, field, ctrl%deltaT)

    case("secondderivatives")
      ! currently only numerical derivatives of forces is implemented

      modeName = "second derivatives"

      ctrl%tDerivs = .true.
      ctrl%tForces = .true.

      call getChildValue(node, "Atoms", buffer2, trim(atomsRange), child=child,&
          & multiple=.true.)
      call getSelectedAtomIndices(child, char(buffer2), geom%speciesNames, geom%species,&
          & ctrl%indDerivAtom)
      if (size(ctrl%indDerivAtom) == 0) then
        call error("No atoms specified for derivatives calculation.")
      end if

      call getChild(node, "MovedAtoms", child, requested=.false.)
      if (associated(child)) then
        if (.not. isContiguousRange(ctrl%indDerivAtom)) then
          call detailedError(child,&
            & "Atoms for calculation of partial Hessian must be a contiguous range.")
        end if
        call getChildValue(child, "", buffer2, child=child2, multiple=.true.)
        call getSelectedAtomIndices(child2, char(buffer2), geom%speciesNames, geom%species, &
           & ctrl%indMovedAtom)
        if (.not. isContiguousRange(ctrl%indMovedAtom)) then
          call detailedError(child2, "MovedAtoms for calculation of partial Hessian must be a &
              & contiguous range.")
        end if
        if (.not. containsAll(ctrl%indDerivAtom, ctrl%indMovedAtom)) then
          call detailedError(child2, "MovedAtoms has indices not contained in Atoms.")
        end if
      else
        ctrl%indMovedAtom = ctrl%indDerivAtom
      end if
      ctrl%nrMoved = size(ctrl%indMovedAtom)

      call getChildValue(node, "Delta", ctrl%deriv2ndDelta, 1.0E-4_dp, &
          & modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, ctrl%deriv2ndDelta)

    case ("velocityverlet")
      ! molecular dynamics

      modeName = "molecular dynamics"

      ctrl%tForces = .true.
      ctrl%tMD = .true.

      call getChildValue(node, "MDRestartFrequency", ctrl%restartFreq, 1)
      call getChildValue(node, "MovedAtoms", buffer2, trim(atomsRange), child=child, &
          &multiple=.true.)
      call getSelectedAtomIndices(child, char(buffer2), geom%speciesNames, geom%species, &
          & ctrl%indMovedAtom)
      ctrl%nrMoved = size(ctrl%indMovedAtom)
      if (ctrl%nrMoved == 0) then
        call error("No atoms specified for molecular dynamics.")
      end if
      call readInitialVelocities(node, ctrl, geom%nAtom)

      call getChildValue(node, "KeepStationary", ctrl%tMDstill,.true.)
      if (ctrl%tMDstill .and. geom%nAtom == 1) then
        call error("Removing translational freedom with only one atom not&
            & possible.")
      end if

      call getChildValue(node, "TimeStep", ctrl%deltaT, modifier=modifier, &
          &child=field)
      call convertUnitHsd(char(modifier), timeUnits, field, ctrl%deltaT)

      call getChildValue(node, "Thermostat", value1, child=child)
      call getNodeName(value1, buffer2)

      thermostat: select case(char(buffer2))
      case ("berendsen")
        ctrl%iThermostat = 2
        ! Read temperature or temperature profiles
        call getChildValue(value1, "Temperature", value2, modifier=modifier, &
            &child=child2)
        call getNodeName(value2, buffer)

        select case(char(buffer))
        case (textNodeName)
          call readTemperature(child2, ctrl)
        case ("temperatureprofile")
          call readTemperatureProfile(value2, char(modifier), ctrl)
        case default
          call detailedError(value2, "Invalid method name.")
        end select

        !call getChildValue(value1, "CouplingStrength", ctrl%wvScale)
        call getChild(value1, "CouplingStrength", child=child2, &
            & requested=.false.)
        if (associated(child2)) then
          call getChildValue(child2, "", ctrl%wvScale)
          call getChild(value1, "Timescale",child=child2,modifier=modifier,&
              &requested=.false.)
          if (associated(child2)) call error("Only Coupling strength OR &
              &Timescale can be set for Berendsen thermostats.")
        else
          call getChild(value1, "Timescale",child=child2,modifier=modifier,&
              &requested=.false.)
          if (associated(child2)) then
            call getChildValue(child2, "", ctrl%wvScale, &
                & modifier=modifier, child=child3)
            call convertUnitHsd(char(modifier), timeUnits, child3, &
                & ctrl%wvScale)
            ctrl%wvScale = ctrl%deltaT / ctrl%wvScale
          else
            call error("Either CouplingStrength or Timescale must be set&
                & for Berendsen thermostats.")
          end if
        end if

      case ("nosehoover")
        ctrl%iThermostat = 3
        ! Read temperature or temperature profiles
        call getChildValue(value1, "Temperature", value2, modifier=modifier, child=child2)
        call getNodeName(value2, buffer)

        select case(char(buffer))
        case (textNodeName)
          call readTemperature(child2, ctrl)
        case ("temperatureprofile")
          call readTemperatureProfile(value2, char(modifier), ctrl)
        case default
          call detailedError(value2, "Invalid method name.")
        end select

        call getChildValue(value1, "CouplingStrength", ctrl%wvScale, modifier=modifier, child=field)
        call convertUnitHsd(char(modifier), freqUnits, field, ctrl%wvScale)

        call getChildValue(value1, "ChainLength", ctrl%nh_npart, 3)
        call getChildValue(value1, "Order", ctrl%nh_nys, 3)
        call getChildValue(value1, "IntegratorSteps", ctrl%nh_nc, 1)

        call getChild(value1, "Restart",  child=child3, requested=.false.)
        if (associated(child3)) then
          allocate(ctrl%xnose(ctrl%nh_npart))
          allocate(ctrl%vnose(ctrl%nh_npart))
          allocate(ctrl%gnose(ctrl%nh_npart))
          call getChildValue(child3,"x",ctrl%xnose)
          call getChildValue(child3,"v",ctrl%vnose)
          call getChildValue(child3,"g",ctrl%gnose)
          ctrl%tInitNHC = .true.
        else
          ctrl%tInitNHC = .false.
        end if

      case ("andersen")
        ctrl%iThermostat = 1
        ! Read temperature or temperature profiles
        call getChildValue(value1, "Temperature", value2, modifier=modifier, child=child2)
        call getNodeName(value2, buffer)

        select case(char(buffer))
        case (textNodeName)
          call readTemperature(child2, ctrl)
        case ("temperatureprofile")
          call readTemperatureProfile(value2, char(modifier), ctrl)
        case default
          call detailedError(value2, "Invalid method name.")
        end select

        call getChildValue(value1, "ReselectProbability", ctrl%wvScale, &
            &child=child3)
        if (ctrl%wvScale <= 0.0_dp .or. ctrl%wvScale > 1.0_dp) then
          call detailedError(child3, &
              &"ReselectProbability must be in the range (0,1]!")
        end if
        call getChildValue(value1, "ReselectIndividually", ctrl%tRescale)

      case ("none")
        ctrl%iThermostat = 0
        allocate(ctrl%tempSteps(1))
        allocate(ctrl%tempValues(1))

        if (ctrl%tReadMDVelocities) then
          ! without a thermostat, if we know the initial velocities, we do not
          ! need a temperature, so just set it to something 'safe'
          ctrl%tempAtom = minTemp
        else
          call readMDInitTemp(value1, ctrl%tempAtom, minTemp)
        end if
      case default
        call getNodeHSDName(value1, buffer2)
        call detailedError(child, "Invalid thermostat '" // char(buffer2) // "'")
      end select thermostat

      if (ctrl%maxRun < -1) then
        call getChildValue(node, "Steps", ctrl%maxRun)
      end if

      call getChildValue(node, "OutputPrefix", buffer2, "geo_end")
      ctrl%outFile = unquote(char(buffer2))

      call getChildValue(node, "Plumed", ctrl%tPlumed, default=.false., child=child)
      if (ctrl%tPlumed .and. .not. withPlumed) then
        call detailedError(child, "Metadynamics can not be used since code has been compiled&
            & without PLUMED support")
      end if

      if (geom%tPeriodic) then

        call getChild(node, "Barostat", child, requested=.false.)
        if (.not. associated(child)) then
          call setChild(node, "Barostat", child)
          ctrl%tBarostat = .false.
          ctrl%pressure = 0.0_dp
          ctrl%BarostatStrength = 0.0_dp
        else
          if (ctrl%nrMoved /= geom%nAtom) then
            call error("Dynamics for a subset of atoms is not currently&
                & possible when using a barostat")
          end if
          call getChildValue(child, "Pressure", ctrl%pressure, &
              & modifier=modifier, child=child2)
          call convertUnitHsd(char(modifier), pressureUnits, child2, &
              & ctrl%pressure)
          call getChild(child, "Coupling", child=child2, requested=.false.)
          if (associated(child2)) then
            call getChildValue(child2, "", ctrl%BarostatStrength)
            call getChild(child, "Timescale",child=child2,modifier=modifier,&
                &requested=.false.)
            if (associated(child2)) call error("Only Coupling strength OR &
                &Timescale can be set for Barostatting.")
          else
            call getChild(child, "Timescale",child=child2,modifier=modifier,&
                &requested=.false.)
            if (associated(child2)) then
              call getChildValue(child2, "", ctrl%BarostatStrength, &
                  & modifier=modifier, child=child3)
              call convertUnitHsd(char(modifier), timeUnits, child3, &
                  & ctrl%BarostatStrength)
              ctrl%BarostatStrength = ctrl%deltaT / ctrl%BarostatStrength
            else
              call error("Either Coupling strength or Timescale must be set&
                  & for Barostatting.")
            end if
          end if
          call getChildValue(child, "Isotropic", ctrl%tIsotropic, .true.)
          ctrl%tBarostat = .true.
        end if
      end if

      if (ctrl%hamiltonian == hamiltonianTypes%dftb) then
        call readXlbomdOptions(node, ctrl%xlbomd)
      end if

      call getInputMasses(node, geom, ctrl%masses)

    case ("socket")
      ! external socket control of the run (once initialised from input)

      modeName = "socket control"

    #:if WITH_SOCKETS
      ctrl%tForces = .true.
      allocate(ctrl%socketInput)
      call getChild(node, 'File', child=child2, requested=.false.)
      call getChild(node, 'Host', child=child3, requested=.false.)
      if (associated(child2) .eqv. associated(child3)) then
        call error('Either Host or File (but not both) must be set for socket&
            & communication')
      end if

      ! File communication
      if (associated(child2)) then
        call getChildValue(child2, "", buffer2)
        ctrl%socketInput%host = unquote(char(buffer2))
        ! zero it to signal to initprogram to use a unix file
        ctrl%socketInput%port = 0
      else
        call getChildValue(child3, "", buffer2)
        ctrl%socketInput%host = unquote(char(buffer2))
        call getChildValue(node, "Port", ctrl%socketInput%port, child=field)
        if (ctrl%socketInput%port <= 0) then
          call detailedError(field, "Invalid port number")
        end if
      end if

      call getChildValue(node, "Protocol", value1, "i-PI", child=child)
      call getNodeName(value1, buffer)
      select case(char(buffer))
      case("i-pi")
        ctrl%socketInput%protocol = IPI_PROTOCOLS%IPI_1
        ! want a file path
        if (ctrl%socketInput%port == 0) then
          call getChildValue(node, "Prefix", buffer2, "/tmp/ipi_")
          sTmp = unquote(char(buffer2))
          ctrl%socketInput%host = trim(sTmp) // trim(ctrl%socketInput%host)
        end if

      case default
        call detailedError(child, "Invalid protocol '" // char(buffer) // "'")
      end select
      call getChildValue(node, "Verbosity", ctrl%socketInput%verbosity, 0)
      call getChildValue(node, "MaxSteps", ctrl%maxRun, 200)

    #:else
      call detailedError(node, "Program had been compiled without socket support")
    #:endif

    case default

      call getNodeHSDName(node, buffer)
      call detailedError(parent, "Invalid driver '" // char(buffer) // "'")

    end select driver

  #:if WITH_TRANSPORT
    if (ctrl%solver%isolver == electronicSolverTypes%OnlyTransport .and. trim(modeName) /= "") then
      call detailederror(node, "transportOnly solver cannot be used with "//trim(modeName))
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
  subroutine commonGeoOptions(node, ctrl, geom, atomsRange, isMaxStepNeeded)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Default range of moving atoms (may be restricted for example by contacts in transport
    !> calculations)
    character(len=*), intent(in) :: atomsRange

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

    call getChildValue(node, "LatticeOpt", ctrl%tLatOpt, .false.)
    if (ctrl%tLatOpt) then
      call getChildValue(node, "Pressure", ctrl%pressure, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), pressureUnits, child, ctrl%pressure)
      call getChildValue(node, "FixAngles", ctrl%tLatOptFixAng, .false.)
      if (ctrl%tLatOptFixAng) then
        call getChildValue(node, "FixLengths", ctrl%tLatOptFixLen, [.false.,.false.,.false.])
      else
        call getChildValue(node, "Isotropic", ctrl%tLatOptIsotropic, .false.)
      end if
      if (isMaxStep) then
        call getChildValue(node, "MaxLatticeStep", ctrl%maxLatDisp, 0.2_dp)
      end if
    end if
    call getChildValue(node, "MovedAtoms", buffer2, trim(atomsRange), child=child, multiple=.true.)
    call getSelectedAtomIndices(child, char(buffer2), geom%speciesNames, geom%species,&
        & ctrl%indMovedAtom)

    ctrl%nrMoved = size(ctrl%indMovedAtom)
    ctrl%tCoordOpt = (ctrl%nrMoved /= 0)
    if (ctrl%tCoordOpt) then
      if (isMaxStep) then
        call getChildValue(node, "MaxAtomStep", ctrl%maxAtomDisp, 0.2_dp)
      end if
    end if
    call getChildValue(node, "MaxForceComponent", ctrl%maxForce, 1e-4_dp, modifier=modifier,&
        & child=field)
    call convertUnitHsd(char(modifier), forceUnits, field, ctrl%maxForce)
    call getChildValue(node, "MaxSteps", ctrl%maxRun, 200)
    call getChildValue(node, "StepSize", ctrl%deltaT, 100.0_dp, modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), timeUnits, field, ctrl%deltaT)
    call getChildValue(node, "OutputPrefix", buffer2, "geo_end")
    ctrl%outFile = unquote(char(buffer2))
    call getChildValue(node, "AppendGeometries", ctrl%tAppendGeo, .false.)
    call readGeoConstraints(node, ctrl, geom%nAtom)
    if (ctrl%tLatOpt) then
      if (ctrl%nrConstr/=0) then
        call error("Lattice optimisation and constraints currently incompatible.")
      end if
      if (ctrl%nrMoved/=0.and.ctrl%nrMoved<geom%nAtom) then
        call error("Subset of optimising atoms not currently possible with lattice optimisation.")
      end if
    end if
    ctrl%isGeoOpt = ctrl%tLatOpt .or. ctrl%tCoordOpt

  end subroutine commonGeoOptions


  !> Extended lagrangian options for XLBOMD
  subroutine readXlbomdOptions(node, input)

    !> node in the input tree
    type(fnode), pointer :: node

    !> extracted settings on exit
    type(TXLBOMDInp), allocatable, intent(out) :: input

    type(fnode), pointer :: pXlbomd, pXlbomdFast, pRoot, pChild
    integer :: nKappa
    logical :: tXlbomdFast

    call getChild(node, 'Xlbomd', pXlbomd, requested=.false.)
    call getChild(node, 'XlbomdFast', pXlbomdFast, requested=.false.)
    if (.not. (associated(pXlbomd) .or. associated(pXlbomdFast))) then
      return
    end if
    if (associated(pXlbomd) .and. associated(pXlbomdFast)) then
      call detailedError(pXlbomdFast, "Blocks 'Xlbomd' and 'XlbomdFast' are&
          & mutually exclusive")
    end if
    if (associated(pXlbomdFast)) then
      tXlbomdFast = .true.
      pRoot => pXlbomdFast
    else
      tXlbomdFast = .false.
      pRoot => pXlbomd
    end if
    allocate(input)
    call getChildValue(pRoot, 'IntegrationSteps', input%nKappa, 5, child=pChild)
    ! Workaround for nagfor 6.1 as the derived type in the array comparison
    ! is not recognized as such
    nKappa = input%nKappa
    if (all([5, 6, 7] /= nKappa)) then
      call detailedError(pChild, 'Invalid number of integration steps (must be&
          & 5, 6 or 7)')
    end if
    call getChildValue(pRoot, 'PreSteps', input%nPreSteps, 0)

    ! Since inverse Jacobian is not enabled, we can set FullSccSteps
    ! to its minimal value (no averaging of inverse Jacobians is done)
    !call getChildValue(child, 'FullSccSteps', input%nFullSccSteps, &
    !    & input%nKappa + 1, child=child2)
    input%nFullSccSteps = input%nKappa + 1
    !if (input%nFullSccSteps < input%nKappa + 1) then
    !  call detailedError(child2, 'Nr. of full SCC steps must be greater by&
    !      & one than integration steps')
    !end if

    if (tXlbomdFast) then
      call getChildValue(pRoot, 'TransientSteps', input%nTransientSteps, 10)
      input%minSccIter = 1
      input%maxSccIter = 1
      ! Dummy value as minSccIter and maxSccIter have been set to 1.
      input%sccTol = 1e-5_dp
      call getChildValue(pRoot, 'Scale', input%scale, 1.0_dp, child=pChild)
      if (input%scale <= 0.0_dp .or. input%scale > 1.0_dp) then
        call detailedError(pChild, 'Scaling value must be in the interval&
            & (0.0, 1.0]')
      end if

      !! Inverse Jacobian is experimental feature so far.
      !call getChildValue(child, "UseJacobianKernel", &
      !    & input%useInverseJacobian, .false.)
      !if (input%useInverseJacobian) then
      !  call getChildValue(child, "ReadJacobianKernel", &
      !      & input%readInverseJacobian, .false.)
      !end if
      input%useInverseJacobian = .false.
      input%readInverseJacobian = .false.
    else
      input%nTransientSteps = 0
      call getChildValue(pRoot, 'MinSccIterations', input%minSCCIter, 1)
      call getChildValue(pRoot, 'MaxSccIterations', input%maxSCCIter, 200)
      if (input%maxSCCIter <= 0) then
        call detailedError(pRoot,"MaxSccIterations must be >= 1");
      end if
      call getChildValue(pRoot, 'SccTolerance', input%sccTol, 1e-5_dp)
      input%scale = 1.0_dp
      input%useInverseJacobian = .false.
      input%readInverseJacobian = .false.
    end if

  end subroutine readXlbomdOptions


  !> Reads geometry constraints
  subroutine readGeoConstraints(node, ctrl, nAtom)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    type(fnode), pointer :: value1, child
    type(string) :: buffer
    type(TListIntR1) :: intBuffer
    type(TListRealR1) :: realBuffer

    call getChildValue(node, "Constraints", value1, "", child=child, allowEmptyValue=.true.)
    call getNodeName2(value1, buffer)
    if (char(buffer) == "") then
      ctrl%nrConstr = 0
    else
      call init(intBuffer)
      call init(realBuffer)
      call getChildValue(child, "", 1, intBuffer, 3, realBuffer)
      ctrl%nrConstr = len(intBuffer)
      allocate(ctrl%conAtom(ctrl%nrConstr))
      allocate(ctrl%conVec(3, ctrl%nrConstr))
      call asVector(intBuffer, ctrl%conAtom)
      if (.not.all(ctrl%conAtom<=nAtom)) then
        call detailedError(node,"Non-existent atom specified in constraint")
      end if
      call asArray(realBuffer, ctrl%conVec)
      call destruct(intBuffer)
      call destruct(realBuffer)
    end if

  end subroutine readGeoConstraints


  !> Reads MD velocities
  subroutine readInitialVelocities(node, ctrl, nAtom)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Total number of all atoms
    integer, intent(in) :: nAtom

    type(fnode), pointer :: value1, child
    type(string) :: buffer, modifier
    type(TListRealR1) :: realBuffer
    integer :: nVelocities
    real(dp), allocatable :: tmpVelocities(:,:)

    call getChildValue(node, "Velocities", value1, "", child=child, &
        & modifier=modifier, allowEmptyValue=.true.)
    call getNodeName2(value1, buffer)
    if (char(buffer) == "") then
      ctrl%tReadMDVelocities = .false.
    else
      call init(realBuffer)
      call getChildValue(child, "", 3, realBuffer, modifier=modifier)
      nVelocities = len(realBuffer)
      if (nVelocities /= nAtom) then
        call detailedError(node, "Incorrect number of specified velocities: " &
            & // i2c(3*nVelocities) // " supplied, " &
            & // i2c(3*nAtom) // " required.")
      end if
      allocate(tmpVelocities(3, nVelocities))
      call asArray(realBuffer, tmpVelocities)
      if (len(modifier) > 0) then
        call convertUnitHsd(char(modifier), VelocityUnits, child, &
            & tmpVelocities)
      end if
      call destruct(realBuffer)
      allocate(ctrl%initialVelocities(3, ctrl%nrMoved))
      ctrl%initialVelocities(:,:) = tmpVelocities(:,ctrl%indMovedAtom(:))
      ctrl%tReadMDVelocities = .true.
    end if

  end subroutine readInitialVelocities


  !> Reads atomic masses from input file, eventually overwriting those in the SK files
  subroutine getInputMasses(node, geo, masses)

    !> relevant node of input data
    type(fnode), pointer :: node

    !> geometry object, which contains atomic species information
    type(TGeometry), intent(in) :: geo

    !> masses to be returned
    real(dp), allocatable, intent(out) :: masses(:)

    type(fnode), pointer :: child, child2, child3, val
    type(fnodeList), pointer :: children
    integer, allocatable :: pTmpI1(:)
    type(string) :: buffer, modifier
    real(dp) :: rTmp
    integer :: ii, jj, iAt

    call getChildValue(node, "Masses", val, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

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
      call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
      call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species, pTmpI1)
      call getChildValue(child2, "MassPerAtom", rTmp, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), massUnits, child, rTmp)
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

    !> Geometry structure
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
  #:else
      call readDFTBHam(node, ctrl, geo, slako, poisson, errStatus)
  #:endif
      @:PROPAGATE_ERROR(errStatus)
    case ("xtb")
  #:if WITH_TRANSPORT
      call readXTBHam(node, ctrl, geo, tp, greendens, poisson, errStatus)
  #:else
      call readXTBHam(node, ctrl, geo, poisson, errStatus)
  #:endif
      @:PROPAGATE_ERROR(errStatus)
    case default
      call detailedError(node, "Invalid Hamiltonian")
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
    type(THybridXcSKTag) :: hybridXcSK

    ctrl%hamiltonian = hamiltonianTypes%dftb

    call readMaxAngularMomentum(node, geo, angShells)

    ! Orbitals and angular momenta for the given shells (once the SK files contain the full
    ! information about the basis, this will be moved to the SK reading routine).
    allocate(slako%orb)
    call setupOrbitals(slako%orb, geo, angShells)

    ! Slater-Koster files
    call getParamSearchPath(searchPath)
    allocate(skFiles(geo%nSpecies, geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      do iSp2 = 1, geo%nSpecies
        call init(skFiles(iSp2, iSp1))
      end do
    end do
    call getChildValue(node, "SlaterKosterFiles", value1, child=child)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case ("type2filenames")
      call getChildValue(value1, "Prefix", buffer2, "")
      prefix = unquote(char(buffer2))
      call getChildValue(value1, "Suffix", buffer2, "")
      suffix = unquote(char(buffer2))
      call getChildValue(value1, "Separator", buffer2, "")
      separator = unquote(char(buffer2))
      call getChildValue(value1, "LowerCaseTypeName", tLower, .false.)
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
            call detailedError(value1, "SK file with generated name '" &
                &// trim(strTmp) // "' does not exist.")
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
          call getChildValue(child, trim(strTmp), lStr, child=child2)
          if (len(lStr) /= len(angShells(iSp1)) * len(angShells(iSp2))) then
            call detailedError(child2, "Incorrect number of Slater-Koster &
                &files")
          end if
          do ii = 1, len(lStr)
            call get(lStr, strTmp, ii)
            call findFile(searchPath, strTmp, strOut)
            if (allocated(strOut)) strTmp = strOut
            inquire(file=strTmp, exist=tExist)
            if (.not. tExist) then
              call detailedError(child2, "SK file '" // trim(strTmp) &
                  &// "' does not exist'")
            end if
            call append(skFiles(iSp2, iSp1), strTmp)
          end do
          call destruct(lStr)
        end do
      end do
    end select

    ! Which repulsive is defined by polynomial? (Default: None)
    allocate(repPoly(geo%nSpecies, geo%nSpecies))
    call getChildValue(node, "PolynomialRepulsive", value1, "", child=child, &
        &list=.true., allowEmptyValue=.true., dummyValue=.true.)
    call getNodeName2(value1, buffer)
    select case (char(buffer))
    case ("")
      repPoly(:,:) = .false.
    case("setforall")
      call getChildValue(value1, "", repPoly(1,1))
      repPoly(:,:) = repPoly(1,1)
    case default
      do iSp1 = 1, geo%nSpecies
        do iSp2 = 1, geo%nSpecies
          strTmp = trim(geo%speciesNames(iSp1)) // "-" &
              &// trim(geo%speciesNames(iSp2))
          call getChildValue(child, trim(strTmp), repPoly(iSp2, iSp1), .false.)
        end do
      end do
      if (.not. all(repPoly .eqv. transpose(repPoly))) then
        call detailedError(value1, "Asymmetric definition (both A-B and B-A must&
            & be defined for using polynomial repulsive)")
      end if
    end select

    call parseChimes(node, ctrl%chimesRepInput)

    ! SCC
    call getChildValue(node, "SCC", ctrl%tSCC, .false.)

    if (ctrl%tSCC) then
      call getChildValue(node, "ShellResolvedSCC", ctrl%tShellResolved, .false.)
    else
      ctrl%tShellResolved = .false.
    end if

    call getChildValue(node, "OldSKInterpolation", ctrl%oldSKInter, .false.)
    if (ctrl%oldSKInter) then
      skInterMeth = skEqGridOld
    else
      skInterMeth = skEqGridNew
    end if

    call parseHybridBlock(node, ctrl%hybridXcInp, geo, skFiles)

    if (.not. allocated(ctrl%hybridXcInp)) then
      call getChild(node, "TruncateSKRange", child, requested=.false.)
      if (associated(child)) then
        call warning("Artificially truncating the SK table, this is normally a bad idea!")
        call SKTruncations(child, rSKCutOff, skInterMeth)
        call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
            & skInterMeth, repPoly, rSKCutOff)
      else
        rSKCutOff = 0.0_dp
        call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
            & skInterMeth, repPoly)
      end if
    else
      call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tShellResolved,&
          & skInterMeth, repPoly, hybridXcSK=hybridXcSK)
      ctrl%hybridXcInp%omega = hybridXcSK%omega
      ctrl%hybridXcInp%camAlpha = hybridXcSK%camAlpha
      ctrl%hybridXcInp%camBeta = hybridXcSK%camBeta
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
      call readSccOptions(node, ctrl, geo)

      ! DFTB hydrogen bond corrections
      call readHCorrection(node, geo, ctrl)

      !> TI-DFTB varibles for Delta DFTB
      call getChild(node, "NonAufbau", child, requested=.false.)
      if (associated(child)) then
        ctrl%isNonAufbau = .true.
        call getChildValue(child, "SpinPurify", ctrl%isSpinPurify, .true.)
        call getChildValue(child, "GroundGuess", ctrl%isGroundGuess, .false.)
        ctrl%tSpin = .true.
        ctrl%t2Component = .false.
        ctrl%nrSpinPol = 0.0_dp
        ctrl%tSpinSharedEf = .false.
      else
        ctrl%isNonAufbau = .false.
      end if

    end if ifSCC

    ! Customize the reference atomic charges for virtual doping
    call readCustomReferenceOcc(node, slako%orb, slako%skOcc, geo, &
        & ctrl%customOccAtoms, ctrl%customOccFillings)

    ! Spin calculation
    if (ctrl%reksInp%reksAlg == reksTypes%noReks  .and. .not.ctrl%isNonAufbau) then
    #:if WITH_TRANSPORT
      call readSpinPolarisation(node, ctrl, geo, tp)
    #:else
      call readSpinPolarisation(node, ctrl, geo)
    #:endif
    end if

    ! temporararily removed until debugged
    !if (.not. ctrl%tscc) then
    !  !! In a non-SCC calculation it is possible to upload charge shifts
    !  !! This is useful if the calculation can jump directly to the Analysis block
    !  call getChildValue(node, "ReadShifts", ctrl%tReadShifts, .false.)
    !end if
    ctrl%tReadShifts = .false.

    ! External fields and potentials
    call readExternal(node, ctrl, geo)

    ! Non-self-consistent spin-orbit coupling
    call readSpinOrbit(node, ctrl, geo, slako%orb)

    ! Electronic solver
  #:if WITH_TRANSPORT
    call readSolver(node, ctrl, geo, tp, greendens, poisson)

    if (tp%taskUpload) then
      ! Initialise variable, but unused
      ctrl%nrChrg =  0.0_dp
    else
      ! Charge
      call getChildValue(node, "Charge", ctrl%nrChrg, 0.0_dp)
    end if
  #:else
    call readSolver(node, ctrl, geo, poisson)

    ! Charge
    call getChildValue(node, "Charge", ctrl%nrChrg, 0.0_dp)
  #:endif

    ! K-Points
    call readKPoints(node, ctrl, geo, ctrl%poorKSampling, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChild(node, "OrbitalPotential", child, requested=.false.)
    if (associated(child)) then
      allocate(ctrl%dftbUInp)
      call getChildValue(child, "Functional", buffer, "fll")
      select case(tolower(char(buffer)))
      case ("fll")
        ctrl%dftbUInp%iFunctional = plusUFunctionals%fll
      case ("psic")
        ctrl%dftbUInp%iFunctional = plusUFunctionals%pSic
      case default
        call detailedError(child,"Unknown orbital functional :"// char(buffer))
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
          call getChildValue(child2,"Shells",li)
          allocate(pTmpI1(len(li)))
          call asArray(li,pTmpI1)
          call append(li1N(iSp1),pTmpI1)
          call append(liN(iSp1),size(pTmpI1))
          deallocate(pTmpI1)
          call destruct(li)
          call getChildValue(child2, "uj", rTmp, 0.0_dp, modifier=modifier, &
              & child=child3)
          call convertUnitHsd(char(modifier), energyUnits, child3, rTmp)
          if (rTmp < 0.0_dp) then
            write(errorStr,"(F12.8)")rTmp
            call detailedError(child2,"Negative value of U-J:"//errorStr)
          end if
          if (rTmp <= 1.0E-10_dp) then
            write(errorStr,"(F12.8)")rTmp
            call detailedError(child2,"Invalid value of U-J, too small: " &
                & //errorStr)
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
    call getChildValue(node, "OnSiteCorrection", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%onSiteElements(slako%orb%mShell, slako%orb%mShell, 2, geo%nSpecies))
      do iSp1 = 1, geo%nSpecies
        call getChildValue(child, trim(geo%speciesNames(iSp1))//"uu",&
            & ctrl%onSiteElements(:slako%orb%nShell(iSp1), :slako%orb%nShell(iSp1), 1, iSp1))
        call getChildValue(child, trim(geo%speciesNames(iSp1))//"ud",&
            & ctrl%onSiteElements(:slako%orb%nShell(iSp1), :slako%orb%nShell(iSp1), 2, iSp1))
      end do
    end if

    ! Dispersion
    call getChildValue(node, "Dispersion", value1, "", child=child, &
        &allowEmptyValue=.true., dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%dispInp)
      call readDispersion(child, geo, ctrl%dispInp, ctrl%nrChrg, ctrl%tSCC)
    end if

    ! Solvation
    call getChildValue(node, "Solvation", value1, "", child=child, &
        &allowEmptyValue=.true., dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%solvInp)
      call readSolvation(child, geo, ctrl%solvInp)
      call getChildValue(value1, "RescaleSolvatedFields", ctrl%isSolvatedFieldRescaled, .true.)
    end if

    ! Electronic constraints
    call getChildValue(node, "ElectronicConstraints", value1, "", child=child,&
        & allowEmptyValue=.true., dummyValue=.true., list=.true.)
    if (associated(value1)) then
      allocate(ctrl%elecConstraintInp)
      call readElecConstraintInput(child, geo, ctrl%elecConstraintInp, ctrl%tSpin)
    end if

    if (ctrl%tLatOpt .and. .not. geo%tPeriodic) then
      call error("Lattice optimisation only applies for periodic structures.")
    end if

  #:if WITH_TRANSPORT
    call readElectrostatics(node, ctrl, geo, tp, poisson)
  #:else
    call readElectrostatics(node, ctrl, geo, poisson)
  #:endif

    ! Third order stuff
    ctrl%t3rd = .false.
    ctrl%t3rdFull = .false.
    if (ctrl%tSCC) then
      call getChildValue(node, "ThirdOrder", ctrl%t3rd, .false.)
      call getChildValue(node, "ThirdOrderFull", ctrl%t3rdFull, .false.)
      if (ctrl%t3rd .and. ctrl%t3rdFull) then
        call detailedError(node, "You must choose either ThirdOrder or&
            & ThirdOrderFull")
      end if
      if (ctrl%t3rd .and. ctrl%tShellResolved) then
        call error("Only full third-order DFTB is compatible with orbital&
            & resolved SCC")
      end if
      if (ctrl%t3rd .or. ctrl%t3rdFull) then
        call getChild(node, 'HubbardDerivs', child, requested=.true.)
        allocate(ctrl%HubDerivs(slako%orb%mShell, geo%nSpecies))
        ctrl%hubDerivs(:,:) = 0.0_dp
        do iSp1 = 1, geo%nSpecies
          nShell = slako%orb%nShell(iSp1)
          if (ctrl%tShellResolved) then
            call getChildValue(child, geo%speciesNames(iSp1),&
                & ctrl%hubDerivs(1:nShell, iSp1))
          else
            call getChildValue(child, geo%speciesNames(iSp1),&
                & ctrl%hubDerivs(1, iSp1))
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
              call getChildValue(node, "HalogenXCorr", ctrl%tHalogenX, .false.)
            end if
          end if
        end block

      end if
    end if

    call readDifferentiation(node, ctrl)

    if (ctrl%tSCC) then
      ! Force type
      call readForceOptions(node, ctrl)
    else
      ctrl%forceType = forceTypes%orig
    end if

    call readCustomisedHubbards(node, geo, slako%orb, ctrl%tShellResolved, ctrl%hubbU)

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

    call getChild(node, "Method", child, requested=.false.)
    if (associated(child)) then
      call getChildValue(child, "", buffer)
      select case(unquote(char(buffer)))
      case default
        call detailedError(child, "Unknown method "//char(buffer)//" for xTB Hamiltonian")
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
      call getChildValue(node, "ParameterFile", value1, "", child=child, &
          &allowEmptyValue=.true., dummyValue=.true.)
      if (associated(value1)) then
        call getChildValue(child, "", buffer)
        paramFile = trim(unquote(char(buffer)))
        call getParamSearchPath(searchPath)
        call findFile(searchPath, paramFile, paramTmp)
        if (allocated(paramTmp)) call move_alloc(paramTmp, paramFile)
        write(stdOut, '(a)') "Using parameter file '"//paramFile//"' for xTB Hamiltonian"
        call ctrl%tbliteInp%setupCalculator(paramFile)
      else
        call detailedError(node, "Either a Method or ParameterFile must be specified for xTB")
      end if
    end if

    call getChildValue(node, "ShellResolvedSCC", ctrl%tShellResolved, .true.)

    ! SCC parameters
    call getChildValue(node, "SCC", ctrl%tSCC, .true.)
    ifSCC: if (ctrl%tSCC) then

      ! get charge mixing options etc.
      call readSccOptions(node, ctrl, geo)

      !> TI-DFTB varibles for Delta DFTB
      call getChild(node, "NonAufbau", child, requested=.false.)
      if (associated(child)) then
        ctrl%isNonAufbau = .true.
        call getChildValue(child, "SpinPurify", ctrl%isSpinPurify, .true.)
        call getChildValue(child, "GroundGuess", ctrl%isGroundGuess, .false.)
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
      call readSpinPolarisation(node, ctrl, geo, tp)
    #:else
      call readSpinPolarisation(node, ctrl, geo)
    #:endif
    end if

    ! temporararily removed until debugged
    !if (.not. ctrl%tscc) then
    !  !! In a non-SCC calculation it is possible to upload charge shifts
    !  !! This is useful if the calculation can jump directly to the Analysis block
    !  call getChildValue(node, "ReadShifts", ctrl%tReadShifts, .false.)
    !end if
    ctrl%tReadShifts = .false.

    ! External fields and potentials
    call readExternal(node, ctrl, geo)

    ! Non-self-consistent spin-orbit coupling
    call ctrl%tbliteInp%setupOrbitals(geo%species, orb)
    call readSpinOrbit(node, ctrl, geo, orb)

    ! Electronic solver
  #:if WITH_TRANSPORT
    call readSolver(node, ctrl, geo, tp, greendens, poisson)

    if (tp%taskUpload) then
      ! Initialise, but unused
      ctrl%nrChrg =  0.0_dp
    else
      ! Charge
      call getChildValue(node, "Charge", ctrl%nrChrg, 0.0_dp)
    end if
  #:else
    call readSolver(node, ctrl, geo, poisson)

    ! Charge
    call getChildValue(node, "Charge", ctrl%nrChrg, 0.0_dp)
  #:endif

    ! K-Points
    call readKPoints(node, ctrl, geo, ctrl%poorKSampling, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    ! Dispersion
    call getChildValue(node, "Dispersion", value1, "", child=child, &
        &allowEmptyValue=.true., dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%dispInp)
      call readDispersion(child, geo, ctrl%dispInp, ctrl%nrChrg, ctrl%tSCC)
    end if

    ! Solvation
    call getChildValue(node, "Solvation", value1, "", child=child, &
        &allowEmptyValue=.true., dummyValue=.true.)
    if (associated(value1)) then
      allocate(ctrl%solvInp)
      call readSolvation(child, geo, ctrl%solvInp)
      call getChildValue(value1, "RescaleSolvatedFields", ctrl%isSolvatedFieldRescaled, .true.)
    end if

    if (ctrl%tLatOpt .and. .not. geo%tPeriodic) then
      call error("Lattice optimisation only applies for periodic structures.")
    end if

  #:if WITH_TRANSPORT
    call readElectrostatics(node, ctrl, geo, tp, poisson)
  #:else
    call readElectrostatics(node, ctrl, geo, poisson)
  #:endif

    ! Third order stuff
    ctrl%t3rd = .true.
    ctrl%t3rdFull = .false.

    call readDifferentiation(node, ctrl)

    if (ctrl%tSCC) then
      ! Force type
      call readForceOptions(node, ctrl)
    else
      ctrl%forceType = forceTypes%orig
    end if

  end subroutine readXTBHam


  !> Reads in settings for spin orbit enabled calculations
  subroutine readSpinOrbit(node, ctrl, geo, orb)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Information about the orbitals of the species/atoms in the system
    class(TOrbitals), intent(in) :: orb

    type(fnode), pointer :: child, child2
    type(string) :: modifier
    integer :: iSp

    call getChild(node, "SpinOrbit", child, requested=.false.)
    if (.not. associated(child)) then
      ctrl%tSpinOrbit = .false.
      allocate(ctrl%xi(0,0))
    else
      if (ctrl%tSpin .and. .not. ctrl%t2Component) then
        call error("Spin-orbit coupling incompatible with collinear spin.")
      end if

      ctrl%tSpinOrbit = .true.
      ctrl%t2Component = .true.

      call getChildValue(child, "Dual", ctrl%tDualSpinOrbit, .true.)

      allocate(ctrl%xi(orb%mShell,geo%nSpecies), source = 0.0_dp)
      do iSp = 1, geo%nSpecies
        call getChildValue(child, geo%speciesNames(iSp), &
            & ctrl%xi(:orb%nShell(iSp),iSp), modifier=modifier, child=child2 )
        call convertUnitHsd(char(modifier), energyUnits, child2,&
            & ctrl%xi(:orb%nShell(iSp),iSp))
      end do
    end if

  end subroutine readSpinOrbit


  !> Read in maximal angular momenta or selected shells
  subroutine readMaxAngularMomentum(node, geo, angShells)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> List containing the angular momenta of the shells
    type(TListIntR1), allocatable, intent(out) :: angShells(:)

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
    call getChild(node, "MaxAngularMomentum", child)
    allocate(angShells(geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      call init(angShells(iSp1))
      call getChildValue(child, geo%speciesNames(iSp1), value1, child=child2)
      call getNodeName(value1, buffer)
      select case(char(buffer))
      case("selectedshells")
        call init(lStr)
        call getChildValue(value1, "", lStr)
        do ii = 1, len(lStr)
          call get(lStr, strTmp, ii)
          strTmp = tolower(unquote(trim(strTmp)))
          if (len_trim(strTmp) > 4 .or. len_trim(strTmp) < 1) then
            call detailedError(value1, "Invalid shell selection '" &
                &// trim(strTmp) &
                &// "'. Nr. of selected shells must be between 1 and 4.")
          end if
          tShellIncl(:) = .false.
          nShell = len_trim(strTmp)
          do jj = 1, nShell
            tmpCh = strTmp(jj:jj)
            tFound = .false.
            do kk = 1, size(shellNames)
              if (tmpCh == trim(shellNames(kk))) then
                if (tShellIncl(kk)) then
                  call detailedError(value1, "Double selection of the same shell&
                      & '" // tmpCh // "' in shell selection block '" &
                      &// trim(strTmp) // "'")
                end if
                tShellIncl(kk) = .true.
                angShell(jj) = kk - 1
                tFound = .true.
                exit
              end if
            end do
            if (.not. tFound) then
              call detailedError(value1, "Invalid shell name '" // tmpCh // "'")
            end if
          end do
          call append(angShells(iSp1), angShell(1:nShell))
        end do
        call destruct(lStr)

      case(textNodeName)
        call getChildValue(child2, "", buffer)
        strTmp = unquote(char(buffer))
        do jj = 1, size(shellNames)
          if (trim(strTmp) == trim(shellNames(jj))) then
            call append(angShells(iSp1), angShellOrdered(:jj))
          end if
        end do
        if (len(angShells(iSp1)) < 1) then
          call detailedError(child2, "Invalid orbital name '" // &
              &trim(strTmp) // "'")
        end if

      case default
        call getNodeHSDName(value1, buffer)
        call detailedError(child2, "Invalid shell specification method '" //&
            & char(buffer) // "'")
      end select
    end do

  end subroutine readMaxAngularMomentum


  !> Setup information about the orbitals of the species/atoms from angShell lists
  subroutine setupOrbitals(orb, geo, angShells)

    !> Information about the orbitals of the species/atoms in the system
    class(TOrbitals), intent(out) :: orb

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> List containing the angular momenta of the shells,
    !> must be inout, since intoArray requires inout arguments
    type(TListIntR1), intent(inout) :: angShells(:)

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
  subroutine readElectrostatics(node, ctrl, geo, tp, poisson)
#:else
  subroutine readElectrostatics(node, ctrl, geo, poisson)
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

    type(fnode), pointer :: value1, child
    type(string) :: buffer

    ctrl%tPoisson = .false.

    ! Read in which kind of electrostatics method to use.
    call getChildValue(node, "Electrostatics", value1, "GammaFunctional", child=child)
    call getNodeName(value1, buffer)

    select case (char(buffer))

    case ("gammafunctional")
    #:if WITH_TRANSPORT
      if (tp%taskUpload .and. ctrl%tSCC) then
        call detailedError(value1, "GammaFunctional not available, if you upload contacts in an SCC&
            & calculation.")
      end if
    #:endif

    case ("poisson")
      if (.not. withPoisson) then
        call detailedError(value1, "Poisson not available as binary was built without the Poisson&
            &-solver")
      end if
      #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
        ctrl%tPoisson = .true.
        #:if WITH_TRANSPORT
          call readPoisson(value1, poisson, geo%tPeriodic, tp, geo%latVecs, ctrl%updateSccAfterDiag)
        #:else
          call readPoisson(value1, poisson, geo%tPeriodic, geo%latVecs, ctrl%updateSccAfterDiag)
        #:endif
      #:endblock

    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Unknown electrostatics '" // char(buffer) // "'")
    end select

  end subroutine readElectrostatics


  !> Spin calculation
#:if WITH_TRANSPORT
  subroutine readSpinPolarisation(node, ctrl, geo, tp)
#:else
  subroutine readSpinPolarisation(node, ctrl, geo)
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

    type(fnode), pointer :: value1, child
    type(string) :: buffer

    call renameChildren(node, "SpinPolarization", "SpinPolarisation")
    call getChildValue(node, "SpinPolarisation", value1, "", child=child, &
        &allowEmptyValue=.true.)
    call getNodeName2(value1, buffer)
    select case(char(buffer))
    case ("")
      ctrl%tSpin = .false.
      ctrl%t2Component = .false.
      ctrl%nrSpinPol = 0.0_dp

    case ("colinear", "collinear")
      ctrl%tSpin = .true.
      ctrl%t2Component = .false.
      call getChildValue(value1, 'UnpairedElectrons', ctrl%nrSpinPol, 0.0_dp)
      call getChildValue(value1, 'RelaxTotalSpin', ctrl%tSpinSharedEf, .false.)
      if (.not. ctrl%tReadChrg) then
        call getInitialSpins(value1, geo, 1, ctrl%initialSpins)
      end if

    case ("noncolinear", "noncollinear")
      ctrl%tSpin = .true.
      ctrl%t2Component = .true.
      if (.not. ctrl%tReadChrg) then
        call getInitialSpins(value1, geo, 3, ctrl%initialSpins)
      end if

    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Invalid spin polarisation type '" //&
          & char(buffer) // "'")
    end select

  end subroutine readSpinPolarisation


  ! External field(s) and potential(s)
  subroutine readExternal(node, ctrl, geo)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

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

    call getChildValue(node, "ElectricField", value1, "", child=child, allowEmptyValue=.true.,&
        & dummyValue=.true., list=.true.)

    ! external applied field
    call getChild(child, "External", child2, requested=.false.)
    if (associated(child2)) then
      allocate(ctrl%electricField)
      ctrl%tMulliken = .true.
      call getChildValue(child2, "Strength", ctrl%electricField%EFieldStrength, modifier=modifier,&
          & child=child3)
      call convertUnitHsd(char(modifier), EFieldUnits, child3, ctrl%electricField%EFieldStrength)
      call getChildValue(child2, "Direction", ctrl%electricField%EfieldVector)
      if (sum(ctrl%electricField%EfieldVector**2) < 1e-8_dp) then
        call detailedError(child2,"Vector too small")
      else
        ctrl%electricField%EfieldVector = ctrl%electricField%EfieldVector&
            & / sqrt(sum(ctrl%electricField%EfieldVector**2))
      end if
      call getChildValue(child2, "Frequency", ctrl%electricField%EFieldOmega, 0.0_dp, &
          & modifier=modifier, child=child3)
      call convertUnitHsd(char(modifier), freqUnits, child3, ctrl%electricField%EFieldOmega)
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
        call getChildValue(child2, "Phase", ctrl%electricField%EfieldPhase, 0)
      end if
    end if

    ctrl%nExtChrg = 0
    if (ctrl%hamiltonian == hamiltonianTypes%dftb) then

      call getChildren(child, "PointCharges", children)
      if (getLength(children) > 0) then
        ! Point charges present
        if (.not.ctrl%tSCC) then
          call error("External charges can only be used in an SCC calculation")
        end if
        call init(lCharges)
        call init(lBlurs)
        ctrl%nExtChrg = 0
        do ii = 1, getLength(children)
          call getItem1(children, ii, child2)
          call getChildValue(child2, "CoordsAndCharges", value1, modifier=modifier, child=child3)
          call getNodeName(value1, buffer)
          select case(char(buffer))
          case (textNodeName)
            call init(lr1)
            call getChildValue(child3, "", 4, lr1, modifier=modifier)
            allocate(tmpR2(4, len(lr1)))
            call asArray(lr1, tmpR2)
            ctrl%nExtChrg = ctrl%nExtChrg + len(lr1)
            call destruct(lr1)
          case ("directread")
            call getChildValue(value1, "Records", ind)
            call getChildValue(value1, "File", buffer2)
            allocate(tmpR2(4, ind))
            call openFile(file, unquote(char(buffer2)), mode="r", iostat=iErr)
            if (iErr /= 0) then
              call detailedError(value1, "Could not open file '"&
                  & // trim(unquote(char(buffer2))) // "' for direct reading" )
            end if
            read(file%unit, *, iostat=iErr) tmpR2
            if (iErr /= 0) then
              call detailedError(value1, "Error during direct reading '"&
                  & // trim(unquote(char(buffer2))) // "'")
            end if
            call closeFile(file)
            ctrl%nExtChrg = ctrl%nExtChrg + ind
          case default
            call detailedError(value1, "Invalid block name")
          end select
          call convertUnitHsd(char(modifier), lengthUnits, child3, tmpR2(1:3,:))
          call append(lCharges, tmpR2)
          call getChildValue(child2, "GaussianBlurWidth", rTmp, 0.0_dp, modifier=modifier,&
              & child=child3)
          if (rTmp < 0.0_dp) then
            call detailedError(child3, "Gaussian blur width may not be negative")
          end if
          call convertUnitHsd(char(modifier), lengthUnits, child3, rTmp)
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
        call detailedError(child, "External charges are not currently supported for this model")
      end if

    end if

    call getChild(node, "AtomSitePotential", child, requested=.false.)
    if (associated(child)) then
      allocate(ctrl%atomicExtPotential)

      call getChild(child, "Net", child2, requested=.false.)
      if (associated(child2)) then
        ! onsites
        ctrl%tNetAtomCharges = .true.
        call init(li)
        call init(lr)
        call getChildValue(child2, "Atoms", li)
        call getChildValue(child2, "Vext", lr, modifier=modifier, child=child3)
        if (len(li) /= len(lr)) then
          call detailedError(child2, "Mismatch in number of sites and potentials")
        end if
        allocate(ctrl%atomicExtPotential%iAtOnSite(len(li)))
        call asArray(li, ctrl%atomicExtPotential%iAtOnSite)
        allocate(ctrl%atomicExtPotential%VextOnSite(len(lr)))
        call asArray(lr, ctrl%atomicExtPotential%VextOnSite)
        call convertUnitHsd(char(modifier), energyUnits, child3, ctrl%atomicExtPotential%VextOnSite)
        call destruct(li)
        call destruct(lr)
      end if

      call getChild(child, "Gross", child2, requested=.false.)
      if (associated(child2)) then
        ! atomic
        call init(li)
        call init(lr)
        call getChildValue(child2, "Atoms", li)
        call getChildValue(child2, "Vext", lr, modifier=modifier, child=child3)
        if (len(li) /= len(lr)) then
          call detailedError(child2, "Mismatch in number of sites and potentials")
        end if
        allocate(ctrl%atomicExtPotential%iAt(len(li)))
        call asArray(li, ctrl%atomicExtPotential%iAt)
        allocate(ctrl%atomicExtPotential%Vext(len(lr)))
        call asArray(lr, ctrl%atomicExtPotential%Vext)
        call convertUnitHsd(char(modifier), energyUnits, child3, ctrl%atomicExtPotential%Vext)
        call destruct(li)
        call destruct(lr)
      end if

      if (.not.allocated(ctrl%atomicExtPotential%iAt)&
          & .and. .not.allocated(ctrl%atomicExtPotential%iAtOnSite)) then
        call detailedError(child, "No atomic potentials specified")
      end if

    end if

  end subroutine readExternal


  !> Filling of electronic levels
  subroutine readFilling(node, ctrl, geo, temperatureDefault)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to test for periodicity
    type(TGeometry), intent(in) :: geo

    !> Default temperature for filling
    real(dp), intent(in) :: temperatureDefault

    type(fnode), pointer :: value1, child, child2, child3, field
    type(string) :: buffer, modifier
    character(lc) :: errorStr

    call getChildValue(node, "Filling", value1, "Fermi", child=child)
    call getNodeName(value1, buffer)

    select case (char(buffer))
    case ("fermi")
      ctrl%iDistribFn = fillingTypes%Fermi ! Fermi function
    case ("methfesselpaxton")
      ! Set the order of the Methfessel-Paxton step function approximation, defaulting to 2nd order
      call getChildValue(value1, "Order", ctrl%iDistribFn, 2)
      if (ctrl%iDistribFn < 1) then
        call getNodeHSDName(value1, buffer)
        write(errorStr, "(A,A,A,I4)")"Unsuported filling mode for '", &
            & char(buffer),"' :",ctrl%iDistribFn
        call detailedError(child, errorStr)
      end if
      ctrl%iDistribFn = fillingTypes%Methfessel + ctrl%iDistribFn
    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Invalid filling method '" //char(buffer)// "'")
    end select

    if (.not. ctrl%tSetFillingTemp) then
      call getChildValue(value1, "Temperature", ctrl%tempElec, temperatureDefault, &
          &modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), energyUnits, field, ctrl%tempElec)
      if (ctrl%tempElec < minTemp) then
        ctrl%tempElec = minTemp
      end if
    end if

    call getChild(value1, "FixedFermiLevel", child=child2, modifier=modifier, requested=.false.)
    ctrl%tFixEf = associated(child2)
    if (ctrl%tFixEf) then
      if (ctrl%tSpin .and. .not.ctrl%t2Component) then
        allocate(ctrl%Ef(2))
      else
        allocate(ctrl%Ef(1))
      end if
      call getChildValue(child2, "", ctrl%Ef, modifier=modifier, child=child3)
      call convertUnitHsd(char(modifier), energyUnits, child3, ctrl%Ef)
    end if

    if (geo%tPeriodic .and. .not.ctrl%tFixEf) then
      call getChildValue(value1, "IndependentKFilling", ctrl%tFillKSep, .false.)
    end if

  end subroutine readFilling


  !> Electronic Solver
#:if WITH_TRANSPORT
  subroutine readSolver(node, ctrl, geo, tp, greendens, poisson)
#:else
  subroutine readSolver(node, ctrl, geo, poisson)
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

    type(fnode), pointer :: value1, child
    type(string) :: buffer, modifier

    integer :: iTmp

    ! Electronic solver
    call getChildValue(node, "Solver", value1, "RelativelyRobust")
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
      call getChildValue(value1, "Sparse", ctrl%solver%elsi%elsiCsr, .false.)
      if (ctrl%solver%elsi%elsiCsr) then
        ctrl%solver%isolver = electronicSolverTypes%elpadm
      else
        ctrl%solver%isolver = electronicSolverTypes%elpa
      end if
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      call getChildValue(value1, "Mode", ctrl%solver%elsi%elpaSolver, 2)
      call getChildValue(value1, "Autotune", ctrl%solver%elsi%elpaAutotune, .false.)
      call getChildValue(value1, "Gpu", ctrl%solver%elsi%elpaGpu, .false., child=child)
      #:if not WITH_GPU
        if (ctrl%solver%elsi%elpaGpu) then
          call detailedError(child, "DFTB+ must be compiled with GPU support in order to enable&
              & the GPU acceleration for the ELPA solver")
        end if
      #:endif

    case ("omm")
      ctrl%solver%isolver = electronicSolverTypes%omm
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      call getChildValue(value1, "nIterationsELPA", ctrl%solver%elsi%ommIterationsElpa, 5)
      call getChildValue(value1, "Tolerance", ctrl%solver%elsi%ommTolerance, 1.0E-10_dp)
      call getChildValue(value1, "Choleskii", ctrl%solver%elsi%ommCholesky, .true.)

    case ("pexsi")
      ctrl%solver%isolver = electronicSolverTypes%pexsi
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
    #:if ELSI_VERSION > 2.5
      call getChildValue(value1, "Method", ctrl%solver%elsi%pexsiMethod, 3)
    #:else
      call getChildValue(value1, "Method", ctrl%solver%elsi%pexsiMethod, 2)
    #:endif
      select case(ctrl%solver%elsi%pexsiMethod)
      case(1)
        iTmp = 60
      case(2)
        iTmp = 20
      case(3)
        iTmp = 30
      end select
      call getChildValue(value1, "Poles", ctrl%solver%elsi%pexsiNPole, iTmp)
      if (ctrl%solver%elsi%pexsiNPole < 10) then
        call detailedError(value1, "Too few PEXSI poles")
      end if
      select case(ctrl%solver%elsi%pexsiMethod)
      case(1)
        if (mod(ctrl%solver%elsi%pexsiNPole,10) /= 0 .or. ctrl%solver%elsi%pexsiNPole > 120) then
          call detailedError(value1, "Unsupported number of PEXSI poles for method 1")
        end if
      case(2,3)
        if (mod(ctrl%solver%elsi%pexsiNPole,5) /= 0 .or. ctrl%solver%elsi%pexsiNPole > 40) then
          call detailedError(value1, "Unsupported number of PEXSI poles for this method")
        end if
      end select
      call getChildValue(value1, "ProcsPerPole", ctrl%solver%elsi%pexsiNpPerPole, 1)
      call getChildValue(value1, "muPoints", ctrl%solver%elsi%pexsiNMu, 2)
      call getChildValue(value1, "SymbolicFactorProcs", ctrl%solver%elsi%pexsiNpSymbo, 1)
      call getChildValue(value1, "SpectralRadius", ctrl%solver%elsi%pexsiDeltaE, 10.0_dp,&
          & modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), energyUnits, child, ctrl%solver%elsi%pexsiDeltaE)

    case ("ntpoly")
      ctrl%solver%isolver = electronicSolverTypes%ntpoly
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      if (ctrl%tSpin) then
        call detailedError(value1, "Solver does not currently support spin polarisation")
      end if
      call getChildValue(value1, "PurificationMethod", ctrl%solver%elsi%ntpolyMethod, 2)
      call getChildValue(value1, "Tolerance", ctrl%solver%elsi%ntpolyTolerance, 1.0E-5_dp)
      call getChildValue(value1, "Truncation", ctrl%solver%elsi%ntpolyTruncation, 1.0E-10_dp)

  #:if WITH_TRANSPORT
    case ("greensfunction")
      ctrl%solver%isolver = electronicSolverTypes%GF
      ! need electronic temperature to be read for this solver:
      call readElectronicFilling(node, ctrl, geo)
      if (tp%defined .and. .not.tp%taskUpload) then
        call detailederror(node, "greensfunction solver cannot be used "// &
            &  "when task = contactHamiltonian")
      end if
      call readGreensFunction(value1, greendens, tp, ctrl%tempElec)
      ! fixEf also avoids checks of total charge in initQFromFile
      ctrl%tFixEf = .true.
    case ("transportonly")
      if (tp%defined .and. .not.tp%taskUpload) then
        call detailederror(node, "transportonly cannot be used when task = contactHamiltonian")
      end if
      call readGreensFunction(value1, greendens, tp, ctrl%tempElec)
      ctrl%solver%isolver = electronicSolverTypes%OnlyTransport
      ctrl%tFixEf = .true.
  #:endif

    case default
      call detailedError(value1, "Unknown electronic solver")

    end select

    if ((ctrl%solver%isolver == electronicSolverTypes%omm .or.&
        & ctrl%solver%isolver == electronicSolverTypes%pexsi ) .and. .not.ctrl%tSpinSharedEf&
        & .and. ctrl%tSpin .and. .not. ctrl%t2Component) then
      call detailedError(value1, "This solver currently requires spin values to be relaxed")
    end if
    if (ctrl%solver%isolver == electronicSolverTypes%pexsi .and. .not.withPEXSI) then
      call error("Not compiled with PEXSI support via ELSI")
    end if
    if (any(ctrl%solver%isolver == [electronicSolverTypes%elpa, electronicSolverTypes%omm,&
        & electronicSolverTypes%pexsi, electronicSolverTypes%ntpoly])) then
      if (.not.withELSI) then
        call error("Not compiled with ELSI supported solvers")
      end if
    end if

    if (any(ctrl%solver%isolver == [electronicSolverTypes%omm, electronicSolverTypes%pexsi,&
        & electronicSolverTypes%ntpoly])) then
      call getChildValue(value1, "Sparse", ctrl%solver%elsi%elsiCsr, .true.)
      if (.not.ctrl%solver%elsi%elsiCsr) then
        if (any(ctrl%solver%isolver == [electronicSolverTypes%pexsi,electronicSolverTypes%ntpoly]))&
            & then
          call getChildValue(value1, "Threshold", ctrl%solver%elsi%elsi_zero_def, 1.0E-15_dp)
        end if
      end if
    end if

  #:if WITH_TRANSPORT
    if (all(ctrl%solver%isolver /= [electronicSolverTypes%GF,electronicSolverTypes%OnlyTransport])&
        & .and. tp%taskUpload) then
      call detailedError(value1, "Eigensolver incompatible with transport calculation&
          & (GreensFunction or TransportOnly required)")
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

    !! Should an additional check be performed if more than one SCC step is requested
    logical :: checkStopHybridCalc

    integer :: ii
    character(lc) :: errorStr

    ! Assume SCC can has usual default number of steps if needed
    poorKSampling = .false.

    ! We can omit any hybrid xc-functional related checks for helical boundary conditions, since
    ! such a calculation will nevertheless be stopped due to the incompatibility of these features
    checkStopHybridCalc = .false.

    ! K-Points
    if (geo%tPeriodic) then
      call getEuclideanKSampling(poorKSampling, checkStopHybridCalc, ctrl, node, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    elseif (geo%tHelical) then
      call getHelicalKSampling(poorKSampling, ctrl, node, geo)
    end if

    call maxSelfConsIterations(node, ctrl, "MaxSCCIterations", ctrl%maxSccIter)
    ! Eventually, perturbation routines should also have restart reads:
    if (ctrl%poorKSampling .and. ctrl%tSCC .and. .not.ctrl%tReadChrg) then
      call warning("It is strongly suggested you use the ReadInitialCharges option.")
    end if

    ! Check if hybrid calculation needs to be stopped due to invalid k-point sampling
    if (checkStopHybridCalc) then
      if (ctrl%maxSccIter == 1) then
        call warning("Restarting a hybrid xc-functional run with what appears to be&
            & a poor k-point sampling that does probably" // NEW_LINE('A') // " not match the&
            & original sampling (however fine for bandstructure calculations).")
      else
        call error("Error while parsing k-point sampling for a hybrid xc-functional&
            & run." // NEW_LINE('A') // "   Only allowed for bandstructure calculations,&
            & i.e. a single SCC iteration.")
      end if
    end if

  end subroutine readKPoints


  !> Set the maximum number of SCC cycles, depending on k-point behaviour
  subroutine maxSelfConsIterations(node, ctrl, label, maxSccIter)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Name of the tag
    character(*), intent(in) :: label

    !> Number of self-consistent iterations
    integer, intent(out) :: maxSccIter

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
      call getChildValue(node, trim(label), maxSccIter, ii)
    end if

    if (ctrl%poorKSampling .and. maxSccIter /= 1) then
      write(warningStr, "(A,I3)") "A self-consistent cycle with these k-points probably will&
          & not correctly calculate many properties, maximum iterations set to:", maxSccIter
      call warning(warningStr)
    end if

  end subroutine maxSelfConsIterations


  !> Tries to infer whether the k-point sampling is restricted to the Gamma-point.
  pure function isGammaOnly(nKPoint, kPoint, kWeight)

    !> Number of k-points for the calculation
    integer, intent(in) :: nKPoint

    !> The k-points for the system
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for the k-points
    real(dp), intent(in) :: kWeight(:)

    !> True, if this appears to be a Gamma-only calculation
    logical :: isGammaOnly

    if (.not. nKPoint == 1) then
      isGammaOnly = .false.
    else
      isGammaOnly = .not. ((.not. all(abs(kPoint(:, 1)) < 1.0e-08_dp))&
          & .or. (.not. abs(kWeight(1)) - 1.0_dp < 1.0e-08_dp))
    end if

  end function isGammaOnly


  !> The k-points in Euclidean space
  subroutine getEuclideanKSampling(poorKSampling, checkStopHybridCalc, ctrl, node, geo, errStatus)

    !> Is this k-point grid usable to integrate properties like the energy, charges, ...?
    logical, intent(out) :: poorKSampling

    !> Should an additional check be performed if more than one SCC step is requested
    logical, intent(out) :: checkStopHybridCalc

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
    real(dp), target :: coeffsAndShifts(3, 4)
    real(dp) :: rTmp3(3)
    type(TListIntR1) :: li1
    type(TListRealR1) :: lr1
    integer, allocatable :: tmpI1(:)
    real(dp), allocatable :: kpts(:,:)

    !! True, if k-points should be reduced by inversion
    logical :: tReduceByInversion

    !! True, if a Gamma-only k-point sampling is requested
    logical :: tGammaOnly

    call getChildValue(node, "KPointsAndWeights", value1, child=child, modifier=modifier)
    call getNodeName(value1, buffer)

    select case(char(buffer))

    case ("supercellfolding")
      poorKSampling = .false.
      if (len(modifier) > 0) then
        call detailedError(child, "No modifier is allowed, if the SupercellFolding scheme is used.")
      end if
      call getChildValue(value1, "", coeffsAndShifts)
      if (abs(determinant33(coeffsAndShifts(:,1:3))) - 1.0_dp < -1e-06_dp) then
        call detailedError(value1, "Determinant of the supercell matrix must be greater than 1")
      end if
      if (any(abs(modulo(coeffsAndShifts(:,1:3) + 0.5_dp, 1.0_dp) - 0.5_dp)&
          & > 1e-06_dp)) then
        call detailedError(value1, "The components of the supercell matrix must be integers.")
      end if
      if (allocated(ctrl%hybridXcInp)) then
        allocate(ctrl%supercellFoldingDiag(3))
        call checkSupercellFoldingMatrix(coeffsAndShifts, errStatus,&
            & supercellFoldingDiagOut=ctrl%supercellFoldingDiag)
        @:PROPAGATE_ERROR(errStatus)
        ctrl%supercellFoldingMatrix = coeffsAndShifts
      end if
      tReduceByInversion = (.not. ctrl%tSpinOrbit)
      call getSuperSampling(coeffsAndShifts(:,1:3), modulo(coeffsAndShifts(:,4), 1.0_dp),&
          & ctrl%kPoint, ctrl%kWeight, reduceByInversion=tReduceByInversion)
      ctrl%nKPoint = size(ctrl%kPoint, dim=2)

    case ("klines")
      ! probably unable to integrate charge for SCC
      poorKSampling = .true.
      call init(li1)
      call init(lr1)
      call getChildValue(value1, "", 1, li1, 3, lr1)
      if (len(li1) < 1) then
        call detailedError(value1, "At least one line must be specified.")
      end if
      allocate(tmpI1(len(li1)))
      allocate(kpts(3, 0:len(lr1)))
      call asVector(li1, tmpI1)
      call asArray(lr1, kpts(:,1:len(lr1)))
      kpts(:,0) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
      call destruct(li1)
      call destruct(lr1)
      if (any(tmpI1 < 0)) then
        call detailedError(value1, "Interval steps must be greater equal to &
            &zero.")
      end if
      ctrl%nKPoint = sum(tmpI1)
      if (ctrl%nKPoint < 1) then
        call detailedError(value1, "Sum of the interval steps must be greater &
            &than zero.")
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
          call detailedError(child, "Invalid modifier: '" // char(modifier) &
              &// "'")
        end select
      end if
      deallocate(tmpI1)
      deallocate(kpts)

    case (textNodeName)

      ! no idea, but assume user knows what they are doing
      poorKSampling = .false.

      call init(lr1)
      call getChildValue(child, "", 4, lr1, modifier=modifier)
      if (len(lr1) < 1) then
        call detailedError(child, "At least one k-point must be defined.")
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
          call detailedError(child, "Invalid modifier: '" // char(modifier) &
              &// "'")
        end select
      end if
      allocate(ctrl%kPoint(3, ctrl%nKPoint))
      allocate(ctrl%kWeight(ctrl%nKPoint))
      ctrl%kPoint(:,:) = kpts(1:3, :)
      ctrl%kWeight(:) = kpts(4, :)
      deallocate(kpts)
    case default
      call detailedError(value1, "Invalid K-point scheme")
    end select

    ! Catch problematic k-point sampling in case this is a hybrid calculation
    checkStopHybridCalc = allocated(ctrl%hybridXcInp) .and. geo%tPeriodic&
        & .and. (char(buffer) /= "supercellfolding") .and. ctrl%tReadChrg

    ! Check for hybrid xc-functional requirements
    tGammaOnly = isGammaOnly(ctrl%nKPoint, ctrl%kPoint, ctrl%kWeight)
    if (.not. tGammaOnly) then
      if (allocated(ctrl%hybridXcInp) .and. geo%tPeriodic&
          & .and. (char(buffer) /= "supercellfolding") .and. (.not. ctrl%tReadChrg)) then
        call detailedError(child, "Error while parsing k-point sampling for a hybrid xc-functional&
            & run. Currently only" // NEW_LINE('A') // "   the supercell folding technique (or any&
            & format specifying the Gamma-point only)" // NEW_LINE('A') // "   is supported.")
      end if
    end if

    ! Hybrid calculations expect the supercell folding coefficients/shifts to be present
    if (allocated(ctrl%hybridXcInp) .and. tGammaOnly&
        & .and. (char(buffer) /= "supercellfolding")) then
      coeffsAndShifts(:,:) = 0.0_dp
      allocate(ctrl%supercellFoldingDiag(3))
      do ii = 1, 3
        coeffsAndShifts(ii, ii) = 1.0_dp
        ctrl%supercellFoldingDiag(ii) = coeffsAndShifts(ii, ii)
      end do
      ctrl%supercellFoldingMatrix = coeffsAndShifts
    end if

  end subroutine getEuclideanKSampling


  !> The k-points for helical boundaries
  subroutine getHelicalKSampling(poorKSampling, ctrl, node, geo)

    !> Is this k-point grid usable to integrate properties like the energy, charges, ...?
    logical, intent(out) :: poorKSampling

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    type(string) :: buffer
    type(fnode), pointer :: value1, child
    type(TListRealR1) :: lr1
    real(dp):: rTmp3(3), rTmp22(2,2)
    integer :: iTmp, iTmp2(2), kk, ii, jj
    real(dp), allocatable :: kPts(:,:)
    character(lc) :: errorStr

    ! assume the user knows what they are doing
    poorKSampling = .false.

    call getChildValue(node, "KPointsAndWeights", value1, child=child)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case ("helicaluniform")
      call getChildValue(value1, "", rTmp3(:2))
      if (abs(modulo(rTmp3(1) + 0.5_dp, 1.0_dp) - 0.5_dp) > 1e-6_dp) then
        call detailedError(value1, "The k-point grid must be integer values.")
      end if
      iTmp = nint(rTmp3(1))
      if (iTmp < 1) then
        call detailedError(node, "Number of grid points must be above 0")
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
        call error("Helical boundaries not yet added for spin-orbit")
      end if
    case ("helicalsampled")
      call getChildValue(value1, "", rTmp22)
      iTmp2 = nint(rTmp22(:,1))
      if (any(abs(iTmp2-rTmp22(:,1)) > 1e-6_dp)) then
        call detailedError(value1, "The k-point grid must be integers.")
      end if
      if (any(iTmp2 < 1)) then
        call detailedError(node, "Number of grid points must be above 0")
      end if
      if (iTmp2(2) > nint(geo%latvecs(3,1))) then
        write(errorStr, '("The k-point grid for the helix rotational operation (",I0,&
            & ") is larger than the rotation order (C_",I0,").")') iTmp2(2), nint(geo%latvecs(3,1))
        call detailedError(node, errorStr)
      end if
      if (mod(nint(geo%latvecs(3,1)),iTmp2(2)) /= 0) then
        write(errorStr, '("The k-point grid for the helix rotational operation (n_k=",I0,&
            & ") is not a divisor of the rotation order (C_",I0,").")') iTmp2(2),&
            & nint(geo%latvecs(3,1))
        call detailedError(node, errorStr)
      end if
      if (abs(rTmp22(2,2) * nint(geo%latvecs(3,1)) - nint(rTmp22(2,2) * nint(geo%latvecs(3,1))))&
          & > epsilon(1.0_dp)) then
        write(errorStr, '("The shift of the k-points along the rotation is incommensurate, it must&
            & be an integer multiple of 1/",I0)') nint(geo%latvecs(3,1))
        call detailedError(node, errorStr)
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
        call error("Helical boundaries not yet added for spin-orbit")
      end if

    case (textNodeName)

      call init(lr1)
      call getChildValue(child, "", 3, lr1)
      if (len(lr1) < 1) then
        call detailedError(child, "At least one k-point must be defined.")
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
        call error("Specified k-value(s) incommensurate with C_n operation.")
      end if
      ! last one is the weight
      ctrl%kWeight(:) = kpts(3, :)
      deallocate(kpts)

    case default
      call detailedError(value1, "Invalid K-point scheme")
    end select

  end subroutine getHelicalKSampling


  !> SCC options that are need for different hamiltonian choices
  subroutine readSccOptions(node, ctrl, geo)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    ctrl%tMulliken = .true.

    call getChildValue(node, "ReadInitialCharges", ctrl%tReadChrg, .false.)
    if (.not. ctrl%tReadChrg) then
      call getInitialCharges(node, geo, ctrl%initialCharges)
    end if

    call getChildValue(node, "SCCTolerance", ctrl%sccTol, 1.0e-5_dp)

    ! temporarily removed until debugged
    ! call getChildValue(node, "WriteShifts", ctrl%tWriteShifts, .false.)
    ctrl%tWriteShifts = .false.

    if (geo%tPeriodic) then
      call getChildValue(node, "EwaldParameter", ctrl%ewaldAlpha, 0.0_dp)
      call getChildValue(node, "EwaldTolerance", ctrl%tolEwald, 1.0e-9_dp)
    end if

    if (geo%tHelical) then
      ! Tolerance for k-points being commensurate with C_n rotation
      call getChildValue(node, "HelicalSymmetryTol", ctrl%helicalSymTol, 1.0E-6_dp)
    end if

    ! self consistency required or not to proceed
    call getChildValue(node, "ConvergentSCCOnly", ctrl%isSccConvRequired, .true.)

  end subroutine readSccOptions


  !> Force evaluation options that are need for different hamiltonian choices
  subroutine readForceOptions(node, ctrl)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    type(fnode), pointer :: child
    type(string) :: buffer

    call getChildValue(node, "ForceEvaluation", buffer, "Traditional", child=child)
    select case (tolower(unquote(char(buffer))))
    case("traditional")
      ctrl%forceType = forceTypes%orig
    case("dynamicst0")
      ctrl%forceType = forceTypes%dynamicT0
    case("dynamics")
      ctrl%forceType = forceTypes%dynamicTFinite
    case default
      call detailedError(child, "Invalid force evaluation method.")
    end select

  end subroutine readForceOptions


  !> Options for truncation of the SK data sets at a fixed distance
  subroutine SKTruncations(node, truncationCutOff, skInterMeth)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> This is the resulting cutoff distance
    real(dp), intent(out) :: truncationCutOff

    !> Method of the sk interpolation
    integer, intent(in) :: skInterMeth

    logical :: tHardCutOff
    type(fnode), pointer :: field
    type(string) :: modifier

    ! Artificially truncate the SK table
    call getChildValue(node, "SKMaxDistance", truncationCutOff, modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), lengthUnits, field, truncationCutOff)

    call getChildValue(node, "HardCutOff", tHardCutOff, .true.)
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
          & goes to 0")
    end if

  end subroutine SKTruncations


  !> Reads initial charges
  subroutine getInitialCharges(node, geo, initCharges)

    !> relevant node in input tree
    type(fnode), pointer :: node

    !> geometry, including atomic type information
    type(TGeometry), intent(in) :: geo

    !> initial atomic charges
    real(dp), allocatable :: initCharges(:)

    type(fnode), pointer :: child, child2, child3, val
    type(fnodeList), pointer :: children
    integer, allocatable :: pTmpI1(:)
    type(string) :: buffer
    real(dp) :: rTmp
    integer :: ii, jj, iAt

    call getChildValue(node, "InitialCharges", val, "", child=child, &
        &allowEmptyValue=.true., dummyValue=.true., list=.true.)

    ! Read either all atom charges, or individual atom specifications
    call getChild(child, "AllAtomCharges", child2, requested=.false.)
    if (associated(child2)) then
      allocate(initCharges(geo%nAtom))
      call getChildValue(child2, "", initCharges)
    else
      call getChildren(child, "AtomCharge", children)
      if (getLength(children) > 0) then
        allocate(initCharges(geo%nAtom))
        initCharges = 0.0_dp
      end if
      do ii = 1, getLength(children)
        call getItem1(children, ii, child2)
        call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
        call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species, pTmpI1)
        call getChildValue(child2, "ChargePerAtom", rTmp)
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
  subroutine getInitialSpins(node, geo, nSpin, initSpins)

    !> relevant node in input data
    type(fnode), pointer :: node

    !> geometry, including atomic information
    type(TGeometry), intent(in) :: geo

    !> number of spin channels
    integer, intent(in) :: nSpin

    !> initial spins on return
    real(dp), allocatable :: initSpins(:,:)

    type(fnode), pointer :: child, child2, child3, val
    type(fnodeList), pointer :: children
    integer, allocatable :: pTmpI1(:)
    type(string) :: buffer
    real(dp), allocatable :: rTmp(:)
    integer :: ii, jj, iAt

    @:ASSERT(nSpin == 1 .or. nSpin == 3)

    call getChildValue(node, "InitialSpins", val, "", child=child, &
        &allowEmptyValue=.true., dummyValue=.true., list=.true.)

    ! Read either all atom spins, or individual spin specifications
    call getChild(child, "AllAtomSpins", child2, requested=.false.)
    if (associated(child2)) then
      allocate(initSpins(nSpin, geo%nAtom))
      call getChildValue(child2, "", initSpins)
    else
      call getChildren(child, "AtomSpin", children)
      if (getLength(children) > 0) then
        allocate(initSpins(nSpin, geo%nAtom))
        initSpins = 0.0_dp
      end if
      allocate(rTmp(nSpin))
      do ii = 1, getLength(children)
        call getItem1(children, ii, child2)
        call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
        call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species, pTmpI1)
        call getChildValue(child2, "SpinPerAtom", rTmp)
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
  subroutine readDifferentiation(node, ctrl)

    !> relevant node in input tree
    type(fnode), pointer, intent(in) :: node

    !> control structure to fill
    type(TControl), intent(inout) :: ctrl


    !> default of a reasonable choice for round off when using a second order finite difference
    !> formula
    real(dp), parameter :: defDelta = epsilon(1.0_dp)**0.25_dp

    type(string) :: buffer, modifier
    type(fnode), pointer :: val, child

    call getChildValue(node, "Differentiation", val, "FiniteDiff",&
        & child=child)
    call getNodeName(val, buffer)
    select case (char(buffer))
    case ("finitediff")
      ctrl%iDerivMethod = diffTypes%finiteDiff
      call getChildValue(val, "Delta", ctrl%deriv1stDelta, defDelta,&
          & modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), lengthUnits, child,&
          & ctrl%deriv1stDelta)
    case ("richardson")
      ctrl%iDerivMethod = diffTypes%richardson
    case default
      call getNodeHSDName(val, buffer)
      call detailedError(child, "Invalid derivative calculation '" &
          & // char(buffer) // "'")
    end select

  end subroutine readDifferentiation


  !> Reads the H corrections (H5, Damp)
  subroutine readHCorrection(node, geo, ctrl)

    !> Node containing the h-bond correction sub-block.
    type(fnode), pointer, intent(in) :: node

    !> Geometry.
    type(TGeometry), intent(in) :: geo

    !> Control structure
    type(TControl), intent(inout) :: ctrl

    type(fnode), pointer :: value1, child, child2
    type(string) :: buffer
    real(dp) :: h5ScalingDef
    integer :: iSp

    ! X-H interaction corrections including H5 and damping
    ctrl%tDampH = .false.
    call getChildValue(node, "HCorrection", value1, "None", child=child)
    call getNodeName(value1, buffer)

    select case (char(buffer))

    case ("none")
      ! nothing to do

    case ("damping")
      ! Switch the correction on
      ctrl%tDampH = .true.
      call getChildValue(value1, "Exponent", ctrl%dampExp)

    case ("h5")
      allocate(ctrl%h5Input)
      associate (h5Input => ctrl%h5Input)
        call getChildValue(value1, "RScaling", h5Input%rScale, 0.714_dp)
        call getChildValue(value1, "WScaling", h5Input%wScale, 0.25_dp)
        allocate(h5Input%elementParams(geo%nSpecies))
        call getChild(value1, "H5Scaling", child2, requested=.false.)
        if (.not. associated(child2)) then
          call setChild(value1, "H5scaling", child2)
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
              & h5ScalingDef)
        end do
        h5Input%speciesNames = geo%speciesNames
      end associate

    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Invalid HCorrection '" // char(buffer) // "'")
    end select

  end subroutine readHCorrection


  !> Reads Slater-Koster files.
  !> Should be replaced with a more sophisticated routine, once the new SK-format has been
  !> established.
  subroutine readSKFiles(skFiles, nSpecies, slako, orb, angShells, orbRes, skInterMeth, repPoly,&
      & truncationCutOff, hybridXcSK, tHyb, tLc, tCam)

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

    !> Distances to artificially truncate tables of SK integrals
    real(dp), intent(in), optional :: truncationCutOff

    !> if calculation range separated then read omega from end of SK file
    type(THybridXcSKTag), intent(inout), optional :: hybridXcSK

    !> True, if global hybrid functional is requested
    logical, intent(in), optional :: tHyb

    !> True, if purely long-range corrected functional is requested
    logical, intent(in), optional :: tLc

    !> True, if CAM range-separation is requested
    logical, intent(in), optional :: tCam

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
            if (.not. present(hybridXcSK)) then
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
                    & hybridXcSK=hybridXcSK)
              elseif (readRep) then
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, iSp1, iSp2,&
                    & splineRepIn=repSplineIn1, hybridXcSK=hybridXcSK)
              else
                call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, hybridXcSK=hybridXcSK)
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
        call checkSKCompElec(skData12, skData21, iSp1, iSp2)
        if (repPoly(iSp1, iSp2)) then
          call checkSKCompRepPoly(repPolyIn1, repPolyIn2, iSp1, iSp2)
        else
          call checkSKCompRepSpline(repSplineIn1, repSplineIn2, iSp1, iSp2)
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
  subroutine checkSKCompElec(skData12, skData21, sp1, sp2)

    !> Slater-Koster integral set for the interaction A-B
    type(TOldSKData), intent(in), target :: skData12(:,:)

    !> Slater-Koster integral set for the interaction B-A
    type(TOldSKData), intent(in), target :: skData21(:,:)

    !> Species number for A (for error messages)
    integer, intent(in) :: sp1

    !> Species number for B (for error messages)
    integer, intent(in) :: sp2

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
          call error(errorStr)
        end if
        if (pSK12%nGrid /= nGrid .or. pSK21%nGrid /= nGrid) then
          write (errorStr, "(A,I2,A,I2)") "Incompatible SK grid lengths for &
              &species pair ", sp1, ", ", sp2
          call error(errorStr)
        end if
      end do
    end do

  end subroutine checkSKCompElec


  !> Checks if the provided repulsive splines for A-B and B-A are compatible
  subroutine checkSKCompRepSpline(repIn1, repIn2, sp1, sp2)

    !> Repulsive spline for interaction A-B
    type(TSplineRepInp), intent(in) :: repIn1

    !> Repulsive spline for interaction B-A
    type(TSplineRepInp), intent(in) :: repIn2

    !> Number of species A (for error messages only)
    integer, intent(in) :: sp1

    !> Number of species B (for error messages only)
    integer, intent(in) :: sp2


    !> Tolerance for the agreement in the repulsive data
    real(dp), parameter :: tolRep = 1.0e-8_dp


    !> string for error return
    character(lc) :: errorStr

    ! Repulsives for A-B and B-A should be the same
    if (size(repIn1%xStart) /= size(repIn2%xStart)) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible nr. of repulsive &
          &intervals for species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if
    if (maxval(abs(repIn1%xStart - repIn2%xStart)) > tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &intervals for species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if
    if (maxval(abs(repIn1%spCoeffs - repIn2%spCoeffs)) > tolRep &
        &.or. maxval(abs(repIn1%spLastCoeffs - repIn2%spLastCoeffs)) &
        &> tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &coefficients for species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if
    if (maxval(abs(repIn1%expCoeffs - repIn2%expCoeffs)) > tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &exp. coefficients for species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if
    if (abs(repIn1%cutoff - repIn2%cutoff) > tolRep) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive spline &
          &cutoffs for species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if

  end subroutine checkSKCompRepSpline


  !> Checks if repulsive polynomials for A-B and B-A are compatible
  subroutine checkSKCompRepPoly(repIn1, repIn2, sp1, sp2)

    !> Repulsive polynomial for interaction A-B
    type(TPolyRepInp), intent(in) :: repIn1

    !> Repulsive polynomial for interaction B-A
    type(TPolyRepInp), intent(in) :: repIn2

    !> Number of species A (for error messages only)
    integer, intent(in) :: sp1

    !> Number of species B (for error messages only)
    integer, intent(in) :: sp2


    !> for error string return
    character(lc) :: errorStr

    if (any(repIn1%polyCoeffs /= repIn2%polyCoeffs)) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive polynomial &
          &coefficients  for the species pair ", sp1, "-", sp2, "."
      call error(errorStr)
    end if
    if (repIn1%cutoff /= repIn2%cutoff) then
      write(errorStr, "(A,I2,A,I2,A)") "Incompatible repulsive cutoffs  &
          &for the species pair ", sp1, "-", sp2, "."
      call error(errorStr)
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
  subroutine readOptions(node, ctrl)

    !> Node to parse
    type(fnode), pointer :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    type(fnode), pointer :: child
    type(string) :: strBuffer
    logical :: tWriteDetailedOutDef

  #:if WITH_SOCKETS
    tWriteDetailedOutDef = .not. allocated(ctrl%socketInput)
  #:else
    tWriteDetailedOutDef = .true.
  #:endif
    call getChildValue(node, "WriteDetailedOut", ctrl%tWriteDetailedOut, tWriteDetailedOutDef)

    call getChildValue(node, "WriteAutotestTag", ctrl%tWriteTagged, .false.)
    call getChildValue(node, "WriteDetailedXML", ctrl%tWriteDetailedXML, .false.)
    call getChildValue(node, "WriteResultsTag", ctrl%tWriteResultsTag, .false.)

    if (.not.(ctrl%tMD.or.ctrl%isGeoOpt.or.allocated(ctrl%geoOpt))) then
      if (ctrl%tSCC) then
        call getChildValue(node, "RestartFrequency", ctrl%restartFreq, 20)
      else
        ctrl%restartFreq = 0
      end if
    end if
    call getChildValue(node, "RandomSeed", ctrl%iSeed, 0, child=child)
    if (ctrl%iSeed < 0) then
      call detailedError(child, "Random seed must be greater or equal zero")
    end if
    call getChildValue(node, "WriteHS", ctrl%tWriteHS, .false.)
    call getChildValue(node, "WriteRealHS", ctrl%tWriteRealHS, .false.)
    call renameChildren(node, "MinimizeMemoryUsage", "MinimiseMemoryUsage")
    call getChildValue(node, "MinimiseMemoryUsage", ctrl%tMinMemory, .false., child=child)
    if (ctrl%tMinMemory) then
      call detailedWarning(child, "Memory minimisation is not working currently, normal calculation&
          & will be used instead")
    end if
    call getChildValue(node, "ShowFoldedCoords", ctrl%tShowFoldedCoord, .false.)
  #:if DEBUG > 0
    call getChildValue(node, "TimingVerbosity", ctrl%timingLevel, -1)
  #:else
    call getChildValue(node, "TimingVerbosity", ctrl%timingLevel, 1)
  #:endif

    if (ctrl%tReadChrg) then
      call getChildValue(node, "ReadChargesAsText", ctrl%tReadChrgAscii, .false.)
    end if

    call getChildValue(node, "WriteCharges", ctrl%tWriteCharges, .true.)
    if (ctrl%tWriteCharges) then
      call getChildValue(node, "WriteChargesAsText", ctrl%tWriteChrgAscii, .false.)
    end if

    ctrl%tSkipChrgChecksum = .false.
    if (.not. ctrl%tFixEf .and. ctrl%tReadChrg) then
      call getChildValue(node, "SkipChargeTest", ctrl%tSkipChrgChecksum, .false.)
    end if

    call readBinaryAccessTypes(node, ctrl%binaryAccessTypes)

  end subroutine readOptions


  !> Reads in dispersion related settings
  subroutine readDispersion(node, geo, input, nrChrg, tSCC)

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

    type(fnode), pointer :: dispModel
    type(string) :: buffer

    call getChildValue(node, "", dispModel)
    call getNodeName(dispModel, buffer)
    select case (char(buffer))
    case ("slaterkirkwood")
      allocate(input%slakirk)
      call readDispSlaKirk(dispModel, geo, input%slakirk)
    case ("lennardjones")
      allocate(input%uff)
      call readDispVdWUFF(dispModel, geo, input%uff)
    case ("dftd3")
      allocate(input%dftd3)
      call readDFTD3(dispModel, geo, input%dftd3)
    case ("simpledftd3")
      allocate(input%sdftd3)
      call readSimpleDFTD3(dispModel, geo, input%sdftd3)
    case ("dftd4")
      allocate(input%dftd4)
      call readDispDFTD4(dispModel, geo, input%dftd4, nrChrg)
    case ("ts")
  #:if WITH_MBD
      allocate(input%mbd)
      call readDispTs(dispModel, input%mbd)
  #:else
      call detailedError(node, "Program must be compiled with the mbd library for TS-dispersion")
  #:endif
    case ("mbd")
  #:if WITH_MBD
      allocate(input%mbd)
      call readDispMbd(dispModel, input%mbd)
  #:else
      call detailedError(node, "Program must be compiled with the mbd library for MBD-dispersion")
  #:endif
    case default
      call detailedError(node, "Invalid dispersion model name.")
    end select

  end subroutine readDispersion


  !> Reads in the dispersion input data for the Slater-Kirkwood dispersion model
  subroutine readDispSlaKirk(node, geo, input)

    !> Node to process
    type(fnode), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the dispersion module on exit
    type(TDispSlaKirkInp), intent(out) :: input

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
    type(TStatus) :: errStatus

    allocate(tmpR2(3, geo%nAtom))
    allocate(input%polar(geo%nAtom))
    allocate(input%rWaals(geo%nAtom))
    allocate(input%charges(geo%nAtom))
    call getChildValue(node, "PolarRadiusCharge", value1, child=child, modifier=modifier)
    call getNodeName(value1, buffer)
    select case (char(buffer))
    case (textNodeName)
      call getChildValue(child, "", tmpR2, modifier=modifier)
      if (len(modifier) > 0) then
        call splitModifier(char(modifier), child, modifiers)
        call convertUnitHsd(char(modifiers(1)), volumeUnits, child, tmpR2(1,:),&
            &.false.)
        call convertUnitHsd(char(modifiers(2)), lengthUnits, child, tmpR2(2,:),&
            &.false.)
        call convertUnitHsd(char(modifiers(3)), chargeUnits, child, tmpR2(3,:),&
            &.false.)
      end if

    case ("hybriddependentpol")
      if (len(modifier) > 0) then
        call detailedError(child, "PolarRadiusCharge is not allowed to carry &
            &a modifier, if the HybridDependentPol method is used.")
      end if
      allocate(rCutoffs(geo%nSpecies))
      allocate(tmp2R2(13, geo%nSpecies))
      do iSp1 = 1, geo%nSpecies
        call getChildValue(value1, geo%speciesNames(iSp1), value2, &
            &child=child2, dummyValue=.true.)
        call getChildValue(child2, "CovalentRadius", rCutoffs(iSp1), &
            &modifier=modifier2, child=child3)
        call convertUnitHsd(char(modifier2), lengthUnits, child3, &
            &rCutoffs(iSp1))
        call renameChildren(child2, "HybridPolarizations", "HybridPolarisations")
        call getChildValue(child2, "HybridPolarisations", tmp2R2(:, iSp1), &
            &modifier=modifier2, child=child3)
        if (len(modifier2) > 0) then
          call splitModifier(char(modifier2), child, modifiers)
          call convertUnitHsd(char(modifiers(1)), volumeUnits, child, &
              &tmp2R2(1:6, iSp1), .false.)
          call convertUnitHsd(char(modifiers(2)), lengthUnits, child, &
              &tmp2R2(7:12, iSp1), .false.)
          call convertUnitHsd(char(modifiers(3)), chargeUnits, child, &
              &tmp2R2(13, iSp1), .false.)
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
      if (errStatus%hasError()) then
        call error(errStatus%message)
      end if
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
      call detailedError(value1, "Invalid method for PolarRadiusCharge.")
    end select

    input%polar(:) = tmpR2(1,:)
    input%rWaals(:) = tmpR2(2,:)
    input%charges(:) = tmpR2(3,:)

  end subroutine readDispSlaKirk


  !> Reads in initialization data for the UFF dispersion model
  subroutine readDispVdWUFF(node, geo, input)

    !> Node to process
    type(fnode), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit
    type(TDispUffInp), intent(out) :: input

    type(string) :: buffer
    type(fnode), pointer :: child, value1, child2
    integer :: iSp
    logical :: found

    call getChildValue(node, "Parameters", value1, child=child)
    allocate(input%distances(geo%nSpecies))
    allocate(input%energies(geo%nSpecies))
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case("uffparameters")
      do iSp = 1, geo%nSpecies
        call getUffValues(geo%speciesNames(iSp), input%distances(iSp), &
            &input%energies(iSp), found)
        if (.not. found) then
          call detailedError(value1, "UFF parameters for species '" // geo&
              &%speciesNames(iSp) // "' not found.")
        end if
      end do
    case default
      call setUnprocessed(value1)
      do iSp = 1, geo%nSpecies
        call getChild(child, geo%speciesNames(iSp), child2)
        call getChildValue(child2, "Distance", input%distances(iSp), &
            &modifier=buffer)
        call convertUnitHsd(char(buffer), lengthUnits, child, &
            &input%distances(iSp))
        call getChildValue(child2, "Energy", input%energies(iSp), &
            &modifier=buffer)
        call convertUnitHsd(char(buffer), energyUnits, child, &
            &input%energies(iSp))
      end do
    end select

  end subroutine readDispVdWUFF


  !> Reads in initialization data for the DFTD3 dispersion module
  subroutine readDFTD3(node, geo, input)

    !> Node to process.
    type(fnode), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit.
    type(TSDFTD3Input), intent(out) :: input

    integer :: iSp
    integer, allocatable :: izpDefault(:)
    type(fnode), pointer :: child, childval
    type(string) :: buffer
    integer, parameter :: d3MaxNum = 94
    logical :: unknownSpecies, threebody

    call getChildValue(node, "Damping", childval, child=child)
    call getNodeName(childval, buffer)
    select case (char(buffer))
    case ("beckejohnson")
      input%dampingFunction = dampingFunction%rational
      call getChildValue(childval, "a1", input%a1)
      call getChildValue(childval, "a2", input%a2)
    case ("zerodamping")
      input%dampingFunction = dampingFunction%zero
      call getChildValue(childval, "sr6", input%sr6)
      call getChildValue(childval, "alpha6", input%alpha6, default=14.0_dp)
    case ("modifiedzerodamping")
      input%dampingFunction = dampingFunction%mzero
      call getChildValue(childval, "sr6", input%sr6)
      call getChildValue(childval, "beta", input%beta)
      call getChildValue(childval, "alpha6", input%alpha6, default=14.0_dp)
    case default
      call getNodeHSDName(childval, buffer)
      call detailedError(child, "Invalid damping method '" // char(buffer) // "'")
    end select
    call getChildValue(node, "s6", input%s6)
    call getChildValue(node, "s8", input%s8)
    call getChildValue(node, "cutoff", input%cutoff, default=sqrt(9000.0_dp), &
        & modifier=buffer, child=child)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoff)
    call getChildValue(node, "cutoffcn", input%cutoffCN, default=40.0_dp, &
        & modifier=buffer, child=child)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoffCN)
    call getChildValue(node, "threebody", threebody, default=.false.)
    input%s9 = merge(1.0_dp, 0.0_dp, threebody)
    ! D3H5 - additional H-H repulsion
    call getChildValue(node, "hhrepulsion", input%hhrepulsion, default=.false.)

    ! Initialize default atomic numbers
    allocate(izpDefault(size(geo%speciesNames)))
    do iSp = 1, size(geo%speciesNames)
      izpDefault(iSp) = symbolToNumber(geo%speciesNames(iSp))
    end do

    ! See if we find user specified overwrites for atomic numbers
    call getChild(node, "AtomicNumbers", child, requested=.false.)
    if (associated(child)) then
      allocate(input%izp(size(geo%speciesNames)))
      call readSpeciesList(child, geo%speciesNames, input%izp, izpDefault)
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
      call detailedError(node, "DFT-D3 does not support all species present")
    end if

  end subroutine readDFTD3


  !> Reads in initialization data for the simple D3 dispersion model.
  subroutine readSimpleDFTD3(node, geo, input)

    !> Node to process.
    type(fnode), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit.
    type(TSimpleDftD3Input), intent(out) :: input

    type(fnode), pointer :: child
    type(string) :: buffer

    call getChildValue(node, "s6", input%s6, default=1.0_dp)
    call getChildValue(node, "s8", input%s8)
    call getChildValue(node, "s10", input%s10, default=0.0_dp)
    call getChildValue(node, "a1", input%a1)
    call getChildValue(node, "a2", input%a2)
    call getChildValue(node, "alpha", input%alpha, default=14.0_dp)
    call getChildValue(node, "weightingFactor", input%weightingFactor, default=4.0_dp)
    call getChildValue(node, "cutoffInter", input%cutoffInter, default=64.0_dp, modifier=buffer,&
        & child=child)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoffInter)

    call readCoordinationNumber(node, input%cnInput, geo, "exp", 0.0_dp)

  end subroutine readSimpleDFTD3


  !> Reads in initialization data for the D4 dispersion model.
  !>
  !> The D4 dispersion model is usually constructed in a failsafe way, so
  !> it only requires to know the damping parameters s8, a1 and a2.
  !> Here we additionally require a s9, since the non-addititive contributions
  !> tend to be expensive especially in the tight-binding context, s9 = 0.0_dp
  !> will disable the calculation.
  subroutine readDispDFTD4(node, geo, input, nrChrg)

    !> Node to process.
    type(fnode), pointer :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Filled input structure on exit.
    type(TDispDftD4Inp), intent(out) :: input

    !> Net charge of the system.
    real(dp), intent(in) :: nrChrg

    integer :: iSp
    integer, allocatable :: izpDefault(:)
    type(fnode), pointer :: value1, child
    type(string) :: buffer
    real(dp), allocatable :: d4Chi(:), d4Gam(:), d4Kcn(:), d4Rad(:)
    integer, parameter :: d4MaxNum = 86
    logical :: unknownSpecies

    call getChildValue(node, "s6", input%s6, default=1.0_dp)
    call getChildValue(node, "s8", input%s8)
    call getChildValue(node, "s9", input%s9)
    call getChildValue(node, "s10", input%s10, default=0.0_dp)
    call getChildValue(node, "a1", input%a1)
    call getChildValue(node, "a2", input%a2)
    call getChildValue(node, "alpha", input%alpha, default=16.0_dp)
    call getChildValue(node, "WeightingFactor", input%weightingFactor, default=6.0_dp)
    call getChildValue(node, "ChargeSteepness", input%chargeSteepness, default=2.0_dp)
    call getChildValue(node, "ChargeScale", input%chargeScale, default=3.0_dp)
    call getChildValue(node, "CutoffInter", input%cutoffInter, default=64.0_dp, modifier=buffer,&
        & child=child)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoffInter)
    call getChildValue(node, "CutoffThree", input%cutoffThree, default=40.0_dp, modifier=buffer,&
        & child=child)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoffThree)

    call getChildValue(node, "ChargeModel", value1, "EEQ", child=child)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(value1, "Unknown method '"//char(buffer)//"' for ChargeModel")
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
      call readEeqModel(value1, input%eeqInput, geo, nrChrg, d4Chi, d4Gam, d4Kcn, d4Rad)
    end select

    ! Initialize default atomic numbers
    allocate(izpDefault(size(geo%speciesNames)))
    do iSp = 1, size(geo%speciesNames)
      izpDefault(iSp) = symbolToNumber(geo%speciesNames(iSp))
    end do

    ! See if we find user specified overwrites for atomic numbers
    call getChild(node, "AtomicNumbers", child, requested=.false.)
    if (associated(child)) then
      allocate(input%izp(size(geo%speciesNames)))
      call readSpeciesList(child, geo%speciesNames, input%izp, izpDefault)
      deallocate(izpDefault)
    else
      call move_alloc(izpDefault, input%izp)
    end if

    call readCoordinationNumber(node, input%cnInput, geo, "Cov", 0.0_dp)

    unknownSpecies = .false.
    do iSp = 1, size(geo%speciesNames)
      if (input%izp(iSp) <= 0 .or. input%izp(iSp) > d4MaxNum) then
        unknownSpecies = .true.
        call warning("Species '"//trim(geo%speciesNames(iSp))// &
          & "' is not supported by DFT-D4")
      end if
    end do
    if (unknownSpecies) then
      call detailedError(node, "DFT-D4 does not support all species present")
    end if

  end subroutine readDispDFTD4


  !> Read settings regarding the EEQ charge model
  subroutine readEeqModel(node, input, geo, nrChrg, kChiDefault, kGamDefault, &
      & kKcnDefault, kRadDefault)

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

    type(fnode), pointer :: value1, child
    type(string) :: buffer

    input%nrChrg = nrChrg

    allocate(input%chi(geo%nSpecies))
    allocate(input%gam(geo%nSpecies))
    allocate(input%kcn(geo%nSpecies))
    allocate(input%rad(geo%nSpecies))

    call getChildValue(node, "Chi", value1, "Defaults", child=child)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' for chi")
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%chi, kChiDefault)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%chi)
    end select

    call getChildValue(node, "Gam", value1, "Defaults", child=child)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' for gam")
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%gam, kGamDefault)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%gam)
    end select

    call getChildValue(node, "Kcn", value1, "Defaults", child=child)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' for kcn")
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%kcn, kKcnDefault)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%kcn)
    end select

    call getChildValue(node, "Rad", value1, "Defaults", child=child)
    call getNodeName(value1, buffer)
    select case(char(buffer))
    case default
      call detailedError(child, "Unknown method '"//char(buffer)//"' for rad")
    case ("defaults")
      call readSpeciesList(value1, geo%speciesNames, input%rad, kRadDefault)
    case ("values")
      call readSpeciesList(value1, geo%speciesNames, input%rad)
    end select

    call getChildValue(node, "Cutoff", input%cutoff, default=40.0_dp, modifier=buffer,&
        & child=child)
    call convertUnitHsd(char(buffer), lengthUnits, child, input%cutoff)

    call getChildValue(node, "EwaldParameter", input%parEwald, 0.0_dp)
    call getChildValue(node, "EwaldTolerance", input%tolEwald, 1.0e-9_dp)

    call readCoordinationNumber(node, input%cnInput, geo, "Erf", 8.0_dp)

  end subroutine readEeqModel


  !> Read in coordination number settings
  subroutine readCoordinationNumber(node, input, geo, cnDefault, cutDefault)

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

    type(fnode), pointer :: value1, value2, child, child2, field
    type(string) :: buffer, modifier
    real(dp), allocatable :: kENDefault(:), kRadDefault(:)

    call getChildValue(node, "CoordinationNumber", value1, cnDefault, child=child)
    call getNodeName(value1, buffer)

    select case(char(buffer))
    case default
      call detailedError(child, "Invalid coordination number type specified")
    case("exp")
      input%cnType = cnType%exp
    case("erf")
      input%cnType = cnType%erf
    case("cov")
      input%cnType = cnType%cov
    case("gfn")
      input%cnType = cnType%gfn
    end select

    call getChildValue(value1, "CutCN", input%maxCN, cutDefault, &
        & child=child2)

    call getChildValue(value1, "Cutoff", input%rCutoff, 40.0_dp, &
        & modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), lengthUnits, field, input%rCutoff)

    allocate(input%en(geo%nSpecies))
    if (input%cnType == cnType%cov) then
      call getChildValue(value1, "Electronegativities", value2, "PaulingEN", child=child2)
      call getNodeName(value2, buffer)
      select case(char(buffer))
      case default
        call detailedError(child2, "Unknown method '" // char(buffer) //&
            & "' to generate electronegativities")
      case("paulingen")
        allocate(kENDefault(geo%nSpecies))
        kENDefault(:) = getElectronegativity(geo%speciesNames)
        call readSpeciesList(value2, geo%speciesNames, input%en, kENDefault)
        deallocate(kENDefault)
      case("values")
        call readSpeciesList(value2, geo%speciesNames, input%en)
      end select
      if (any(input%en <= 0.0_dp)) then
        call detailedError(value1, "Electronegativities are not defined for all species")
      end if
    else
      ! array is not used, but should still be populated with dummies
      input%en(:) = 0.0_dp
    end if

    allocate(input%covRad(geo%nSpecies))
    call getChildValue(value1, "Radii", value2, "CovalentRadiiD3", child=child2)
    call getNodeName(value2, buffer)
    select case(char(buffer))
    case default
      call detailedError(child2, "Unknown method '"//char(buffer)//"' to generate radii")
    case("covalentradiid3")
      allocate(kRadDefault(geo%nSpecies))
      kRadDefault(:) = getCovalentRadius(geo%speciesNames)
      call readSpeciesList(value2, geo%speciesNames, input%covRad, kRadDefault)
      deallocate(kRadDefault)
    case("values")
      call readSpeciesList(value2, geo%speciesNames, input%covRad)
    end select

    if (any(input%covRad <= 0.0_dp)) then
      call detailedError(value1, "Covalent radii are not defined for all species")
    end if

  end subroutine readCoordinationNumber


#:if WITH_MBD

  !> Reads in settings for Tkatchenko-Scheffler dispersion
  subroutine readDispTs(node, input)

    !> data to parse
    type(fnode), pointer, intent(in) :: node

    !> control data coming back
    type(TDispMbdInp), intent(out) :: input

    type(string) :: buffer
    type(fnode), pointer :: child

    input%method = 'ts'
    call getChild(node, "EnergyAccuracy", child, requested=.false.)
    if (associated(child)) then
      call detailedWarning(child, "The energy accuracy setting will be ignored as it is not&
          & supported/need by libMBD any more")
    end if
    call getChild(node, "ForceAccuracy", child, requested=.false.)
    if (associated(child)) then
      call detailedWarning(child, "The force accuracy setting will be ignored as it is not&
          & supported/need by libMBD any more")
    end if
    call getChildValue(node, "Damping", input%ts_d, default=(input%ts_d))
    call getChildValue(node, "RangeSeparation", input%ts_sr, default=(input%ts_sr))
    call getChildValue(node, "ReferenceSet", buffer, 'ts', child=child)
    input%vdw_params_kind = tolower(unquote(char(buffer)))
    call checkManyBodyDispRefName(input%vdw_params_kind, child)
    call getChildValue(node, "LogLevel", input%log_level, default=(input%log_level))
  end subroutine readDispTs


  !> Reads in many-body dispersion settings
  subroutine readDispMbd(node, input)

    !> data to parse
    type(fnode), pointer, intent(in) :: node

    !> control data coming back
    type(TDispMbdInp), intent(out) :: input

    type(string) :: buffer
    type(fnode), pointer :: child

    input%method = 'mbd-rsscs'
    call getChildValue(node, "Beta", input%mbd_beta, input%mbd_beta)
    call getChildValue(node, "NOmegaGrid", input%n_omega_grid, default=(input%n_omega_grid))
    call getChildValue(node, "KGrid", input%k_grid)
    call getChildValue(node, "KGridShift", input%k_grid_shift, default=(input%k_grid_shift))
    call getChildValue(node, "ReferenceSet", buffer, 'ts', child=child)
    input%vdw_params_kind = tolower(unquote(char(buffer)))
    call checkManyBodyDispRefName(input%vdw_params_kind, child)
    call getChildValue(node, "LogLevel", input%log_level, default=(input%log_level))

  end subroutine readDispMbd


  !> Check the dispersion label matches allowed cases
  subroutine checkManyBodyDispRefName(name, node)

    !> Label
    character(*), intent(in) :: name

    !> data tree for error usage
    type(fnode), pointer, intent(in) :: node

    if (name /= 'ts' .and. name /= 'tssurf') then
      call detailedError(node, 'Invalid reference set name for TS/MBD-dispersion')
    end if

  end subroutine checkManyBodyDispRefName

#:endif

  !> reads in value of temperature for MD with sanity checking of the input
  subroutine readTemperature(node, ctrl)

    !> data to parse
    type(fnode), pointer :: node

    !> control data coming back
    type(TControl), intent(inout) :: ctrl

    type(string) :: modifier

    allocate(ctrl%tempSteps(1))
    allocate(ctrl%tempValues(1))
    allocate(ctrl%tempMethods(1))
    ctrl%tempMethods(1) = 1
    ctrl%tempSteps(1) = 1
    call getChildValue(node, "", ctrl%tempValues(1), modifier=modifier)
    call convertUnitHsd(char(modifier), energyUnits, node, ctrl%tempValues(1))
    if (ctrl%tempValues(1) < 0.0_dp) then
      call detailedError(node, "Negative temperature.")
    end if
    if (ctrl%tempValues(1) < minTemp) then
      ctrl%tempValues(1) = minTemp
    end if

  end subroutine readTemperature


  !> reads a temperature profile for MD with sanity checking of the input
  subroutine readTemperatureProfile(node, modifier, ctrl)

    !> parser node containing the relevant part of the user input
    type(fnode), pointer :: node

    !> unit modifier for the profile
    character(len=*), intent(in) :: modifier

    !> Control structure to populate
    type(TControl), intent(inout) :: ctrl

    type(TListString) :: ls
    type(TListIntR1) :: li1
    type(TListRealR1) :: lr1
    character(len=20), allocatable :: tmpC1(:)
    integer :: ii
    logical :: success

    call init(ls)
    call init(li1)
    call init(lr1)
    call getChildValue(node, "", ls, 1, li1, 1, lr1)
    if (len(ls) < 1) then
      call detailedError(node, "At least one annealing step must be &
          &specified.")
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
      call detailedError(node, "Invalid annealing method name '" // trim(tmpC1(ii)) // "'.")
    end do

    if (any(ctrl%tempSteps < 0)) then
      call detailedError(node, "Step values must not be negative.")
    end if

    ii = sum(ctrl%tempSteps)
    if (ii < 1) then
      call detailedError(node, "Sum of steps in the profile must be &
          &greater than zero.")
    end if
    ctrl%maxRun = ii - 1

    if (any(ctrl%tempValues < 0.0_dp)) then
      call detailedError(node, "Negative temperature.")
    end if

    call convertUnitHsd(modifier, energyUnits, node, ctrl%tempValues)
    if (any(ctrl%tempValues < minTemp)) then
      ctrl%tempValues = max(ctrl%tempValues, minTemp)
    end if
    deallocate(tmpC1)

  end subroutine readTemperatureProfile


  !> Reads the excited state data block
  subroutine readExcited(node, geo, ctrl)

    !> Node to parse
    type(fnode), pointer :: node

    !> geometry object, which contains atomic species information
    type(TGeometry), intent(in) :: geo

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    type(fnode), pointer :: child
    type(fnode), pointer :: child2, child3
    type(fnode), pointer :: value
    type(string) :: buffer, modifier

    ! Linear response stuff
    call getChild(node, "Casida", child, requested=.false.)

    if (associated(child)) then

      allocate(ctrl%lrespini)
      ctrl%lrespini%tPrintEigVecs = .false.

      if (ctrl%tSpin) then
        ctrl%lrespini%sym = ' '
      else
        call getChildValue(child, "Symmetry", buffer, child=child2)
        select case (unquote(char(buffer)))
        case ("Singlet" , "singlet")
          ctrl%lrespini%sym = 'S'
        case ("Triplet" , "triplet")
          ctrl%lrespini%sym = 'T'
        case ("Both" , "both")
          ctrl%lrespini%sym = 'B'
        case default
          call detailedError(child2, "Invalid symmetry value '"  // char(buffer) // &
              & "' (must be 'Singlet', 'Triplet' or 'Both').")
        end select
      end if

      call getChildValue(child, "NrOfExcitations", ctrl%lrespini%nexc)

      call getChild(child, "StateOfInterest", child2, requested=.false.)
      if (.not. associated(child2)) then
        ctrl%lrespini%nstat = 0
        call setChildValue(child, "StateOfInterest", 0)
      else
        call getChildValue(child2, "", buffer)
        if (tolower(unquote(char(buffer))) == "brightest") then
          if (ctrl%lrespini%sym /= "S" .or. ctrl%tSpin) then
            call detailedError(child2, "Brightest mode only allowed for spin unpolarised singlet&
                & excitations.")
          end if
          ctrl%lrespini%nstat = -1
        else
          call getChildValue(child2, "", ctrl%lrespini%nstat)
          if (ctrl%lrespini%nstat > ctrl%lrespini%nexc) then
            call detailedError(child2, "Invalid value, must be within range of NrOfExcitations")
          elseif (ctrl%lrespini%sym == "B" .and. ctrl%lrespini%nstat /= 0) then
            call detailedError(child2, "You cannot specify a particular excited state if symmetry&
                & is 'B'")
          end if
        end if
      end if

      call getChildValue(child, "EnergyWindow", ctrl%lrespini%energyWindow, 0.0_dp, &
          & modifier=modifier, child=child2)
      ctrl%lrespini%tEnergyWindow = ctrl%lrespini%energyWindow /= 0.0_dp
      call convertUnitHsd(char(modifier), energyUnits, child2, ctrl%lrespini%energyWindow)
      call getChildValue(child, "OscillatorWindow", ctrl%lrespini%oscillatorWindow, 0.0_dp, &
          & modifier=modifier,  child=child2)
      ctrl%lrespini%tOscillatorWindow = ctrl%lrespini%oscillatorWindow /= 0.0_dp
      call convertUnitHsd(char(modifier), dipoleUnits, child2, ctrl%lrespini%oscillatorWindow)
      call getChildValue(child, "CacheCharges", ctrl%lrespini%tCacheCharges, default=.true.)
      call getChildValue(child, "WriteMulliken", ctrl%lrespini%tMulliken, default=.false.)
      call getChildValue(child, "WriteCoefficients", ctrl%lrespini%tCoeffs, default=.false.)
      ctrl%lrespini%tGrndState = .false.
      if (ctrl%lrespini%tCoeffs) then
        call getChildValue(child, "TotalStateCoeffs", ctrl%lrespini%tGrndState, .false.)
      end if
      call getChildValue(child, "WriteEigenvectors", ctrl%lrespini%tPrintEigVecs, .false.)
      call getChildValue(child, "WriteDensityMatrix", ctrl%lrespini%tWriteDensityMatrix, .false.)
      call getChildValue(child, "WriteXplusY", ctrl%lrespini%tXplusY, default=.false.)
      call getChildValue(child, "StateCouplings", ctrl%lrespini%indNACouplings, default=[0, 0])
      if (all(ctrl%lrespini%indNACouplings == 0)) then
        ctrl%lrespini%tNaCoupling = .false.
      else
        ctrl%lrespini%tNaCoupling = .true.
      end if
      call getChildValue(child, "WriteSPTransitions", ctrl%lrespini%tSPTrans, default=.false.)
      call getChildValue(child, "WriteTransitions", ctrl%lrespini%tTrans, default=.false.)
      call getChildValue(child, "WriteTransitionDipole", ctrl%lrespini%tTradip, default=.false.)
      if (allocated(ctrl%hybridXcInp)) then
        call getChildValue(child, "WriteTransitionCharges", ctrl%lrespini%tTransQ, default=.false.)
      end if
      ctrl%lrespini%iLinRespSolver = linRespSolverTypes%None

      call renameChildren(child, "Diagonalizer", "Diagonaliser")
      call getChildValue(child, "Diagonaliser", child2)
      call getNodeName(child2, buffer)
      select case(char(buffer))
      case ("arpack")
        if (.not. withArpack) then
          call detailedError(child2, 'This DFTB+ binary has been compiled without support for&
              & linear response calculations using the ARPACK/ngARPACK library.')
        end if
        call getChildValue(child2, "WriteStatusArnoldi", ctrl%lrespini%tArnoldi, default=.false.)
        call getChildValue(child2, "TestArnoldi", ctrl%lrespini%tDiagnoseArnoldi, default=.false.)
        ctrl%lrespini%iLinRespSolver = linRespSolverTypes%Arpack
      case ("stratmann")
        ctrl%lrespini%iLinRespSolver = linRespSolverTypes%Stratmann
        call getChildValue(child2, "SubSpaceFactor", ctrl%lrespini%subSpaceFactorStratmann, 20)
      case default
        call detailedError(child2, "Invalid diagonaliser method '" // char(buffer) // "'")
      end select

      call getChildValue(child, "OptimiserCI", child2, "", child=child3, allowEmptyValue=.true.)
      if (associated(child2)) then
        call getNodeName(child2, buffer)
        select case(char(buffer))
        case ("bearpark")
          ctrl%lrespini%isCIopt = .true.
          call getChildValue(child2, "EnergyShift", ctrl%lrespini%energyShiftCI,&
              & modifier=modifier, default=0.0_dp)
          call convertUnitHsd(char(modifier), energyUnits, child, ctrl%lrespini%energyShiftCI)
        case default
          call detailedError(child2, "Invalid optimiser method '" // char(buffer) // "'")
        end select
      else
        ctrl%lrespini%isCIopt = .false.
      end if

      if (ctrl%tForces .or. ctrl%tPrintForces) then
        call getChildValue(child, "ExcitedStateForces", ctrl%tCasidaForces, default=.true.)
      end if

    end if

    !pp-RPA
    call getChild(node, "PP-RPA", child, requested=.false.)

    if (associated(child)) then

      allocate(ctrl%pprpa)

      if (ctrl%tSpin) then
        ctrl%pprpa%sym = ' '
      else
        call getChildValue(child, "Symmetry", buffer, child=child2)
        select case (unquote(char(buffer)))
        case ("Singlet" , "singlet")
          ctrl%pprpa%sym = 'S'
        case ("Triplet" , "triplet")
          ctrl%pprpa%sym = 'T'
        case ("Both" , "both")
          ctrl%pprpa%sym = 'B'
        case default
          call detailedError(child2, "Invalid symmetry value '"  // char(buffer) // &
              & "' (must be 'Singlet', 'Triplet' or 'Both').")
        end select
      end if

      call getChildValue(child, "NrOfExcitations", ctrl%pprpa%nexc)

      call getChildValue(child, "HHubbard", value, child=child2)
      allocate(ctrl%pprpa%hhubbard(geo%nSpecies))
      call readSpeciesList(child2, geo%speciesNames, ctrl%pprpa%hhubbard)

      call getChildValue(child, "TammDancoff", ctrl%pprpa%tTDA, default=.false.)

      call getChild(child, "NrOfVirtualStates", child2, requested=.false.)
      if (.not. associated(child2)) then
        ctrl%pprpa%nvirtual = 0
        ctrl%pprpa%tConstVir = .false.
        call setChildValue(child, "NrOfVirtualStates", 0)
      else
        call getChildValue(child2, "", ctrl%pprpa%nvirtual)
        ctrl%pprpa%tConstVir = .true.
      end if

    end if

  end subroutine readExcited


  !> Reads the analysis block
#:if WITH_TRANSPORT
  subroutine readAnalysis(node, ctrl, geo, orb, transpar, tundos)
#:else
  subroutine readAnalysis(node, ctrl, geo)
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

      call getChildValue(node, "ProjectStates", val, "", child=child, &
          & allowEmptyValue=.true., list=.true.)
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
          call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
          call getSelectedAtomIndices(child3, char(buffer), geo%speciesNames, geo%species, pTmpI1)
          call append(ctrl%iAtInRegion, pTmpI1)
          call getChildValue(child2, "ShellResolved", &
              & ctrl%tShellResInRegion(iReg), .false., child=child3)
          if (ctrl%tShellResInRegion(iReg)) then
            if (.not. all(geo%species(pTmpI1) == geo%species(pTmpI1(1)))) then
              call detailedError(child3, "Shell resolved PDOS only allowed for &
                  &regions where all atoms belong to the same species")
            end if
          end if
          call getChildValue(child2, "OrbitalResolved", &
              & ctrl%tOrbResInRegion(iReg), .false., child=child3)
          if (ctrl%tOrbResInRegion(iReg)) then
            if (.not. all(geo%species(pTmpI1) == geo%species(pTmpI1(1)))) then
              call detailedError(child3, "Orbital resolved PDOS only allowed for &
                  &regions where all atoms belong to the same species")
            end if
          end if
          deallocate(pTmpI1)
          write(strTmp, "('region',I0)") iReg
          call getChildValue(child2, "Label", buffer, trim(strTmp))
          ctrl%RegionLabel(iReg) = unquote(char(buffer))
        end do
      end if
      call destroyNodeList(children)

      call renameChildren(node, "Localize", "Localise")
      call getChild(node, "Localise", child=val, requested=.false.)
      if (associated(val)) then
        ctrl%tLocalise = .true.
        call getChild(val, "PipekMezey", child=child2, requested=.false.)
        if (associated(child2)) then
          allocate(ctrl%pipekMezeyInp)
          associate(inp => ctrl%pipekMezeyInp)
            call getChildValue(child2, "MaxIterations", inp%maxIter, 100)
            tPipekDense = .true.
            if (.not. geo%tPeriodic) then
              call getChildValue(child2, "Dense", tPipekDense, .false.)
              if (.not. tPipekDense) then
                call init(lr1)
                call getChild(child2, "SparseTolerances", child=child3, requested=.false.)
                if (associated(child3)) then
                  call getChildValue(child3, "", 1, lr1)
                  if (len(lr1) < 1) then
                    call detailedError(child2, "Missing values of tolerances.")
                  end if
                  allocate(inp%sparseTols(len(lr1)))
                  call asVector(lr1, inp%sparseTols)
                else
                  allocate(inp%sparseTols(4))
                  inp%sparseTols = [0.1_dp, 0.01_dp, 1.0E-6_dp, 1.0E-12_dp]
                  call setChildValue(child2, "SparseTolerances", inp%sparseTols)
                end if
                call destruct(lr1)
              end if
            end if
            if (tPipekDense) then
              call getChildValue(child2, "Tolerance", inp%tolerance, 1.0E-4_dp)
            end if
          end associate
        else
          call detailedError(val, "No localisation method chosen")
        end if
      end if

      call getChildValue(node, "WriteEigenvectors", ctrl%tPrintEigVecs, .false.)

    #:if WITH_SOCKETS
      tWriteBandDatDef = .not. allocated(ctrl%socketInput)
    #:else
      tWriteBandDatDef = .true.
    #:endif

      call getChildValue(node, "WriteBandOut", ctrl%tWriteBandDat, tWriteBandDatDef)

      call getChild(node, "Polarisability", child=child, requested=.false.)
      call getChild(node, "ResponseKernel", child=child2, requested=.false.)
      if (associated(child) .or. associated(child2)) then
        allocate(ctrl%perturbInp)
      end if

      ! electric field polarisability of system
      call getChild(node, "Polarisability", child=child, requested=.false.)
      if (associated(child)) then
        ctrl%perturbInp%isEPerturb = .true.
        call freqRanges(child, ctrl%perturbInp%dynEFreq)
      end if

      call getChild(node, "ResponseKernel", child=child, requested=.false.)
      if (associated(child)) then
        ctrl%perturbInp%isRespKernelPert = .true.
        if (ctrl%tSCC) then
          call getChildValue(child, "RPA", ctrl%perturbInp%isRespKernelRPA, .false.)
        else
          ctrl%perturbInp%isRespKernelRPA = .true.
        end if
        call freqRanges(child, ctrl%perturbInp%dynKernelFreq)
      end if

      if (allocated(ctrl%perturbInp)) then
        call getChildValue(node, "PertubDegenTol", ctrl%perturbInp%tolDegenDFTBPT, 128.0_dp,&
            & child=child)
        if (ctrl%perturbInp%tolDegenDFTBPT < 1.0_dp) then
          call detailedError(child, "Perturbation degeneracy tolerance must be above 1x")
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
          call getChildValue(node, "PerturbEta", ctrl%perturbInp%etaFreq, 1.0E-8_dp, child=child)
          if (ctrl%perturbInp%etaFreq < epsilon(0.0_dp)) then
            call detailedError(child, "Imaginary constant for finite frequency perturbation too&
                & small")
          end if
        end if
      end if

      if (allocated(ctrl%perturbInp)) then
        call maxSelfConsIterations(node, ctrl, "MaxPerturbIter", ctrl%perturbInp%maxPerturbIter)
        if (ctrl%tScc) then
          call getChildValue(node, "PerturbSccTol", ctrl%perturbInp%perturbSccTol, 1.0e-5_dp)
          ! self consistency required, or not, to proceed with perturbation
          call getChildValue(node, "ConvergedPerturb", ctrl%perturbInp%isPerturbConvRequired,&
              & .true.)
        end if
      end if

    end if

    if (tHaveDensityMatrix) then

      ! Is this compatible with Poisson solver use?
      call readElectrostaticPotential(node, geo, ctrl)

      call getChildValue(node, "MullikenAnalysis", ctrl%tPrintMulliken, .true.)
      if (ctrl%tPrintMulliken) then
        call getChildValue(node, "WriteNetCharges", ctrl%tPrintNetAtomCharges, default=.false.)
        if (ctrl%tPrintNetAtomCharges) then
          ctrl%tNetAtomCharges = .true.
        end if
        call getChild(node, "CM5", child, requested=.false.)
        if (associated(child)) then
          allocate(ctrl%cm5Input)
          call readCM5(child, ctrl%cm5Input, geo)
        end if
      end if
      call getChildValue(node, "AtomResolvedEnergies", ctrl%tAtomicEnergy, .false.)

      if (allocated(ctrl%solvInp)) then
        call getChildValue(node, "writeCosmoFile", ctrl%tWriteCosmoFile, &
            & allocated(ctrl%solvInp%cosmoInp), child=child)
        if (ctrl%tWriteCosmoFile .and. .not.allocated(ctrl%solvInp%cosmoInp)) then
          call detailedError(child, "Cosmo file can only be written for Cosmo calculations")
        end if
      end if

      call getChildValue(node, "PrintForces", ctrl%tPrintForces, .false.)

    else

      ctrl%tPrintMulliken = .false.
      ctrl%tAtomicEnergy = .false.
      ctrl%tPrintForces = .false.

    end if


  #:if WITH_TRANSPORT
    call getChild(node, "TunnelingAndDOS", child, requested=.false.)
    if (associated(child)) then
      if (.not.transpar%defined) then
        call error("Block TunnelingAndDos requires Transport block.")
      end if
      if (.not.transpar%taskUpload) then
        call error("Block TunnelingAndDos not compatible with task=contactHamiltonian")
      end if
      if (.not. allocated(orb)) then
        call error("Orbital information from SK-files missing (xTB Hamiltonian not compatible&
            & with transport yet)")
      end if
      call readTunAndDos(child, orb, geo, tundos, transpar, ctrl%tempElec)
    else
      if (ctrl%solver%isolver == electronicSolverTypes%OnlyTransport) then
        call detailedError(node, "The TransportOnly solver requires a TunnelingAndDos block to be&
            & present.")
      end if
    endif
  #:endif

  end subroutine readAnalysis


  !> Frequency ranges for response calculations
  subroutine freqRanges(node, frequencies)

    !> Node to parse
    type(fnode), pointer :: node

    !> Frequencies, 0 being static
    real(dp), allocatable, intent(inout) :: frequencies(:)

    type(TListReal) :: lr
    type(fnode), pointer :: child, child2
    type(string) :: modifier
    integer :: nFreq, iFreq, jFreq
    real(dp) :: tmp3R(3)
    logical :: isStatic

    call getChildValue(node, "Static", isStatic, .true.)
    if (isStatic) then
      call growFreqArray(frequencies, 1)
      ! should already be zero, but just in case:
      frequencies(:) = 0.0_dp
    end if

    call getChild(node, "Frequencies", child=child, modifier=modifier, requested=.false.)
    if (associated(child)) then
      call init(lr)
      call getChildValue(child, "", lr, child=child2, modifier=modifier)
      nFreq = len(lr)
      if (nFreq > 0) then
        if (allocated(frequencies)) then
          iFreq = size(frequencies)
        else
          iFreq = 0
        end if
        call growFreqArray(frequencies, nFreq)
        call asArray(lr, frequencies(iFreq+1:iFreq+nFreq))
        call convertUnitHsd(char(modifier),freqUnits, child, frequencies(iFreq+1:iFreq+nFreq))
      end if
      call destruct(lr)
      if (any(frequencies < 0.0_dp)) then
        call detailedError(child2, "Negative driving frequency requested")
      end if
    end if

    call getChild(node, "FrequencyRange", child=child, modifier=modifier, requested=.false.)
    if (associated(child)) then
      call init(lr)
      call getChildValue(child, "", lr, child=child2, modifier=modifier)
      if (len(lr) == 3) then
        call asArray(lr, tmp3R)
        call convertUnitHsd(char(modifier), freqUnits, child, tmp3R)
        if (any(tmp3R(:2) < 0.0_dp)) then
          call detailedError(child, "Negative values in dynamic frequency range.")
        end if
        if (abs(tmp3R(3)) <= epsilon(0.0_dp)) then
          call detailedError(child, "Increase step size in dynamic frequency range.")
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
        call detailedError(child,"Malformed frequency range.")
      end if
      call destruct(lr)
      if (any(frequencies < 0.0_dp)) then
        call detailedError(child2, "Negative driving frequency requested")
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
  subroutine readLaterAnalysis(node, ctrl)

    !> Node to parse
    type(fnode), pointer :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl


    logical :: tPrintEigVecs

    tPrintEigVecs = ctrl%tPrintEigVecs
    if (allocated(ctrl%lrespini)) tPrintEigvecs = tPrintEigvecs .or. ctrl%lrespini%tPrintEigVecs
    if (tPrintEigVecs) then
      call getChildValue(node, "EigenvectorsAsText", ctrl%tPrintEigVecsTxt, .false.)
    end if

  end subroutine readLaterAnalysis


  !> Read in hamiltonian settings that are influenced by those read from REKS{}, electronDynamics{}
  subroutine readLaterHamiltonian(hamNode, ctrl, driverNode, geo)

    !> Hamiltonian node to parse
    type(fnode), pointer :: hamNode

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> Geometry driver node to parse
    type(fnode), pointer :: driverNode

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    type(fnode), pointer :: value1, value2, child, child2
    type(string) :: buffer, buffer2
    type(TListRealR1) :: lr1

    if (ctrl%reksInp%reksAlg == reksTypes%noReks) then

      if (ctrl%tSCC) then

        call getChildValue(hamNode, "Mixer", value1, "Broyden", child=child)
        call getNodeName(value1, buffer)
        select case(char(buffer))

        case ("broyden")

          ctrl%iMixSwitch = mixerTypes%broyden
          call getChildValue(value1, "MixingParameter", ctrl%almix, 0.2_dp)
          call getChildValue(value1, "InverseJacobiWeight", ctrl%broydenOmega0, 0.01_dp)
          call getChildValue(value1, "MinimalWeight", ctrl%broydenMinWeight, 1.0_dp)
          call getChildValue(value1, "MaximalWeight", ctrl%broydenMaxWeight, 1.0e5_dp)
          call getChildValue(value1, "WeightFactor", ctrl%broydenWeightFac, 1.0e-2_dp)

        case ("anderson")

          ctrl%iMixSwitch = mixerTypes%anderson
          call getChildValue(value1, "MixingParameter", ctrl%almix, 0.05_dp)
          call getChildValue(value1, "Generations", ctrl%iGenerations, 4)
          call getChildValue(value1, "InitMixingParameter", ctrl%andersonInitMixing, 0.01_dp)
          call getChildValue(value1, "DynMixingParameters", value2, "", child=child,&
              & allowEmptyValue=.true.)
          call getNodeName2(value2, buffer2)
          if (char(buffer2) == "") then
            ctrl%andersonNrDynMix = 0
          else
            call init(lr1)
            call getChildValue(child, "", 2, lr1, child=child2)
            if (len(lr1) < 1) then
              call detailedError(child2, "At least one dynamic mixing parameter must be defined.")
            end if
            ctrl%andersonNrDynMix = len(lr1)
            allocate(ctrl%andersonDynMixParams(2, ctrl%andersonNrDynMix))
            call asArray(lr1, ctrl%andersonDynMixParams)
            call destruct(lr1)
          end if
          call getChildValue(value1, "DiagonalRescaling", ctrl%andersonOmega0, 1.0e-2_dp)

        case ("simple")

          ctrl%iMixSwitch = mixerTypes%simple
          call getChildValue(value1, "MixingParameter", ctrl%almix, 0.05_dp)

        case("diis")

          ctrl%iMixSwitch = mixerTypes%diis
          call getChildValue(value1, "InitMixingParameter", ctrl%almix, 0.2_dp)
          call getChildValue(value1, "Generations", ctrl%iGenerations, 6)
          call getChildValue(value1, "UseFromStart", ctrl%tFromStart, .true.)

        case default

          call getNodeHSDName(value1, buffer)
          call detailedError(child, "Invalid mixer '" // char(buffer) // "'")

        end select

      end if

      if (ctrl%tMD) then
        if (ctrl%iThermostat /= 0) then
          call getChildValue(driverNode, "Thermostat", child, child=child2)
          if (ctrl%reksInp%reksAlg == reksTypes%noReks) then
            call getChildValue(child, "AdaptFillingTemp", ctrl%tSetFillingTemp, .false.)
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
        call readElectronicFilling(hamNode, ctrl, geo)
      end if

    end if hamNeedsT

  end subroutine readLaterHamiltonian


  !> Parses for electronic filling temperature (should only read if not either REKS or electron
  !> dynamics from a supplied density matrix)
  subroutine readElectronicFilling(hamNode, ctrl, geo)

    !> Relevant node in input tree
    type(fnode), pointer :: hamNode

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure to test for periodicity
    type(TGeometry), intent(in) :: geo

    select case(ctrl%hamiltonian)
    case(hamiltonianTypes%xtb)
      call readFilling(hamNode, ctrl, geo, 300.0_dp*Boltzmann)
    case(hamiltonianTypes%dftb)
      call readFilling(hamNode, ctrl, geo, 0.0_dp)
    end select

  end subroutine readElectronicFilling


  !> Reads W values if required by settings in the Hamiltonian or the excited state
  subroutine readSpinConstants(hamNode, geo, orb, ctrl)

    !> node for Hamiltonian data
    type(fnode), pointer :: hamNode

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> control structure
    type(TControl), intent(inout) :: ctrl

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

      call getChild(hamNode, "SpinConstants", child)
      if (ctrl%hamiltonian == hamiltonianTypes%xtb) then
        call getChildValue(child, "ShellResolvedSpin", tShellResolvedW, .true.)
      else
        if (.not.ctrl%tShellResolved) then
          call getChildValue(child, "ShellResolvedSpin", tShellResolvedW, .false.)
        else
          tShellResolvedW = .true.
        end if
      end if

      if (tShellResolvedW) then
        ! potentially unique values for each shell
        do iSp1 = 1, geo%nSpecies
          call getChildValue(child, geo%speciesNames(iSp1),&
              & ctrl%spinW(:orb%nShell(iSp1), :orb%nShell(iSp1), iSp1))
        end do
      else
        ! only one value per atom
        do iSp1 = 1, geo%nSpecies
          call getChildValue(child, geo%speciesNames(iSp1),ctrl%spinW(1, 1, iSp1))
          ctrl%spinW(:orb%nShell(iSp1), :orb%nShell(iSp1), iSp1) =&
              & ctrl%spinW(1, 1, iSp1)
        end do
      end if
    end if

  end subroutine readSpinConstants


  !> Reads customised Hubbard U values that over-ride the SK file values
  subroutine readCustomisedHubbards(node, geo, orb, tShellResolvedScc, hubbU)

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

    type(fnode), pointer :: child, child2
    integer :: iSp1

    call renameChildren(node, "CustomizedHubbards", "CustomisedHubbards")
    call getChild(node, "CustomisedHubbards", child, requested=.false.)
    if (associated(child)) then
      allocate(hubbU(orb%mShell, geo%nSpecies))
      hubbU(:,:) = 0.0_dp
      do iSp1 = 1, geo%nSpecies
        call getChild(child, geo%speciesNames(iSp1), child2, requested=.false.)
        if (.not. associated(child2)) then
          cycle
        end if
        if (tShellResolvedScc) then
          call getChildValue(child2, "", hubbU(:orb%nShell(iSp1), iSp1))
        else
          call getChildValue(child2, "", hubbU(1, iSp1))
          hubbU(:orb%nShell(iSp1), iSp1) = hubbU(1, iSp1)
        end if
      end do
    end if

  end subroutine readCustomisedHubbards


  !> Reads the electron dynamics block
  subroutine readElecDynamics(node, input, geom, masses)

    !> input data to parse
    type(fnode), pointer :: node

    !> ElecDynamicsInp instance
    type(TElecDynamicsInp), intent(inout) :: input

    !> geometry of the system
    type(TGeometry), intent(in) :: geom

    !> masses to be returned
    real(dp), allocatable, intent(inout) :: masses(:)

    type(fnode), pointer :: value1, child
    type(string) :: buffer, buffer2, modifier
    logical :: ppRangeInvalid, tNeedFieldStrength
    real (dp) :: defPpRange(2)
    logical :: defaultWrite

    call getChildValue(node, "Steps", input%steps)
    call getChildValue(node, "TimeStep", input%dt, modifier=modifier, child=child)
    call convertUnitHsd(char(modifier), timeUnits, child, input%dt)

    call getChildValue(node, "Populations", input%tPopulations, .false.)
    call getChildValue(node, "WriteFrequency", input%writeFreq, 50)
    call getChildValue(node, "Restart", input%tReadRestart, .false.)
    if (input%tReadRestart) then
      call getChildValue(node, "RestartFromAscii", input%tReadRestartAscii, .false.)
    end if
    call getChildValue(node, "WriteRestart", input%tWriteRestart, .true.)
    if (input%tWriteRestart) then
      call getChildValue(node, "WriteAsciiRestart", input%tWriteRestartAscii, .false.)
    end if
    call getChildValue(node, "RestartFrequency", input%restartFreq, max(input%Steps / 10, 1))
    call getChildValue(node, "Forces", input%tForces, .false.)
    call getChildValue(node, "WriteBondEnergy", input%tBondE, .false.)
    call getChildValue(node, "WriteBondPopulation", input%tBondP, .false.)
    call getChildValue(node, "WriteAtomicEnergies", input%tWriteAtomEnergies, .false.)
    call getChildValue(node, "Pump", input%tPump, .false.)
    call getChildValue(node, "FillingsFromFile", input%tFillingsFromFile, .false.)

    if (input%tPump) then
      call getChildValue(node, "PumpProbeFrames", input%tdPPFrames)
      defPpRange = [0.0_dp, input%steps * input%dt]
      call getChildValue(node, "PumpProbeRange", input%tdPpRange, defPprange, modifier=modifier,&
          & child=child)
      call convertUnitHsd(char(modifier), timeUnits, child, input%tdPpRange)

      ppRangeInvalid = (input%tdPpRange(2) <= input%tdPpRange(1))&
          & .or. (input%tdPprange(1) < defPpRange(1))&
          & .or. (input%tdPpRange(2) > defPpRange(2))
      if (ppRangeInvalid) then
        call detailederror(child, "Wrong definition of PumpProbeRange, either incorrect order&
            & or outside of simulation time range")
      end if
    end if

    call getChildValue(node, "Probe", input%tProbe, .false.)
    if (input%tPump .and. input%tProbe) then
      call detailedError(child, "Pump and probe cannot be simultaneously true.")
    end if

    call getChildValue(node, "EulerFrequency", input%eulerFreq, 0)
    if ((input%eulerFreq < 50) .and. (input%eulerFreq > 0)) then
      call detailedError(child, "Wrong number of Euler steps, should be above 50")
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
    call getChildValue(node, "Perturbation", value1, "None", child=child)
    call getNodeName(value1, buffer)
    select case(char(buffer))

    case ("kick")
      input%pertType = pertTypes%kick
      call renameChildren(value1, "PolarizationDirection", "PolarisationDirection")
      call getChildValue(value1, "PolarisationDirection", buffer2)
      input%polDir = directionConversion(unquote(char(buffer2)), value1)

      call getChildValue(value1, "SpinType", buffer2, "Singlet")
      select case(unquote(char(buffer2)))
      case ("singlet", "Singlet")
        input%spType = tdSpinTypes%singlet
      case ("triplet", "Triplet")
        input%spType = tdSpinTypes%triplet
      case default
        call detailedError(value1, "Unknown spectrum spin type " // char(buffer2))
      end select

      defaultWrite = .false.

    case ("laser")
      input%pertType = pertTypes%laser
      call renameChildren(value1, "PolarizationDirection", "PolarisationDirection")
      call getChildValue(value1, "PolarisationDirection", input%reFieldPolVec)
      call renameChildren(value1, "ImagPolarizationDirection", "ImagPolarisationDirection")
      call getChildValue(value1, "ImagPolarisationDirection", input%imFieldPolVec, &
          & [0.0_dp, 0.0_dp, 0.0_dp])
      call getChildValue(value1, "LaserEnergy", input%omega, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), energyUnits, child, input%omega)
      call getChildValue(value1, "Phase", input%phase, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), angularUnits, child, input%phase)
      call getChildValue(value1, "ExcitedAtoms", buffer, "1:-1", child=child, multiple=.true.)
      call getSelectedAtomIndices(child, char(buffer), geom%speciesNames, geom%species,&
          & input%indExcitedAtom)

      input%nExcitedAtom = size(input%indExcitedAtom)
      if (input%nExcitedAtom == 0) then
        call error("No atoms specified for laser excitation.")
      end if

      defaultWrite = .true.

    case ("kickandlaser")
      input%pertType = pertTypes%kickAndLaser
      call getChildValue(value1, "KickPolDir", buffer2)
      input%polDir = directionConversion(unquote(char(buffer2)), value1)
      call getChildValue(value1, "SpinType", input%spType, tdSpinTypes%singlet)
      call getChildValue(value1, "LaserPolDir", input%reFieldPolVec)
      call getChildValue(value1, "LaserImagPolDir", input%imFieldPolVec, [0.0_dp, 0.0_dp, 0.0_dp])
      call getChildValue(value1, "LaserEnergy", input%omega, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), energyUnits, child, input%omega)
      call getChildValue(value1, "Phase", input%phase, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), angularUnits, child, input%phase)
      call getChildValue(value1, "LaserStrength", input%tdLaserField, modifier=modifier,&
          & child=child)
      call convertUnitHsd(char(modifier), EFieldUnits, child, input%tdLaserField)

      call getChildValue(value1, "ExcitedAtoms", buffer, "1:-1", child=child, multiple=.true.)
      call getSelectedAtomIndices(child, char(buffer), geom%speciesNames, geom%species,&
          & input%indExcitedAtom)
      input%nExcitedAtom = size(input%indExcitedAtom)
      if (input%nExcitedAtom == 0) then
        call error("No atoms specified for laser excitation.")
      end if

      defaultWrite = .false.

    case ("none")
      input%pertType = pertTypes%noTDPert
      tNeedFieldStrength = .false.

      defaultWrite = .true.

    case default
      call detailedError(child, "Unknown perturbation type " // char(buffer))
    end select

    if (tNeedFieldStrength) then
      call getChildValue(node, "FieldStrength", input%tdfield, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), EFieldUnits, child, input%tdfield)
    end if

    call getChildValue(node, "WriteEnergyAndCharges", input%tdWriteExtras, defaultWrite)

    !! Different envelope functions
    call getChildValue(node, "EnvelopeShape", value1, "Constant")
    call getNodeName(value1, buffer)
    select case(char(buffer))

    case("constant")
      input%envType = envTypes%constant

    case("gaussian")
      input%envType = envTypes%gaussian
      call getChildValue(value1, "Time0", input%time0, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), timeUnits, child, input%Time0)

      call getChildValue(value1, "Time1", input%time1, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), timeUnits, child, input%Time1)

    case("sin2")
      input%envType = envTypes%sin2
      call getChildValue(value1, "Time0", input%time0, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), timeUnits, child, input%Time0)

      call getChildValue(value1, "Time1", input%time1, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), timeUnits, child, input%Time1)

    case("fromfile")
      input%envType = envTypes%fromFile
      call getChildValue(value1, "Time0", input%time0, 0.0_dp, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), timeUnits, child, input%Time0)

    case default
      call detailedError(value1, "Unknown envelope shape " // char(buffer))
    end select

    !! Non-adiabatic molecular dynamics
    call getChildValue(node, "IonDynamics", input%tIons, .false.)
    if (input%tIons) then
      call getChildValue(node, "MovedAtoms", buffer, "1:-1", child=child, multiple=.true.)
      call getSelectedAtomIndices(child, char(buffer), geom%speciesNames, geom%species,&
          & input%indMovedAtom)

      input%nMovedAtom = size(input%indMovedAtom)
      call readInitialVelocitiesNAMD(node, input, geom%nAtom)
      if (input%tReadMDVelocities) then
        ! without a thermostat, if we know the initial velocities, we do not need a temperature, so
        ! just set it to something 'safe'
        input%tempAtom = minTemp
      else
        if (.not. input%tReadRestart) then
          ! previously lower limit was minTemp:
          call readMDInitTemp(node, input%tempAtom, 0.0_dp)
        end if
        call getInputMasses(node, geom, masses)
      end if
    end if

  end subroutine readElecDynamics


  !> Read in initial ion temperature for simple MD
  subroutine readMDInitTemp(node, tempAtom, minimumTemp)

    !> input data to parse
    type(fnode), pointer :: node

    !> Ionic temperature
    real(dp), intent(out) :: tempAtom

    !> Lowest possible ion temperature
    real(dp), intent(in) :: minimumTemp

    type(fnode), pointer :: child
    type(string) :: modifier

    call getChildValue(node, "InitialTemperature", tempAtom, modifier=modifier, child=child)
    if (tempAtom < 0.0_dp) then
      call detailedError(node, "Negative temperature")
    end if
    call convertUnitHsd(char(modifier), energyUnits, node, tempAtom)
    tempAtom = max(tempAtom, minimumTemp)

  end subroutine readMDInitTemp


  !> Converts direction label text string into corresponding numerical value
  function directionConversion(direction, node) result(iX)

    !> Direction label
    character(*), intent(in) :: direction

    !> input tree for error return
    type(fnode), pointer :: node

    !> direction indicator (1 - 4) for (x,y,z,all)
    integer :: iX

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
          & // ". Must be x, y, z or all.")
    end select

  end function directionConversion


  !> Reads MD velocities
  subroutine readInitialVelocitiesNAMD(node, input, nAtom)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> ElecDynamicsInp object structure to be filled
    type(TElecDynamicsInp), intent(inout) :: input

    !> Total number of all atoms
    integer, intent(in) :: nAtom

    type(fnode), pointer :: value1, child
    type(string) :: buffer, modifier
    type(TListRealR1) :: realBuffer
    integer :: nVelocities
    real(dp), pointer :: tmpVelocities(:,:)

    call getChildValue(node, "Velocities", value1, "", child=child, &
         & modifier=modifier, allowEmptyValue=.true.)
    call getNodeName2(value1, buffer)
    if (char(buffer) == "") then
       input%tReadMDVelocities = .false.
    else
       call init(realBuffer)
       call getChildValue(child, "", 3, realBuffer, modifier=modifier)
       nVelocities = len(realBuffer)
       if (nVelocities /= nAtom) then
          call detailedError(node, "Incorrect number of specified velocities: " &
               & // i2c(3*nVelocities) // " supplied, " &
               & // i2c(3*nAtom) // " required.")
       end if
       allocate(tmpVelocities(3, nVelocities))
       call asArray(realBuffer, tmpVelocities)
       if (len(modifier) > 0) then
          call convertUnitHsd(char(modifier), VelocityUnits, child, &
               & tmpVelocities)
       end if
       call destruct(realBuffer)
       allocate(input%initialVelocities(3, input%nMovedAtom))
       input%initialVelocities(:,:) = tmpVelocities(:, input%indMovedAtom(:))
       input%tReadMDVelocities = .true.
    end if

  end subroutine readInitialVelocitiesNAMD


#:if WITH_TRANSPORT
  !> Read geometry information for transport calculation
  subroutine readTransportGeometry(root, geom, transpar)

    !> Root node containing the current block
    type(fnode), pointer :: root

    !> geometry of the system, which may be modified for some types of calculation
    type(TGeometry), intent(inout) :: geom

    !> Parameters of the transport calculation
    type(TTransPar), intent(inout) :: transpar

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
    call getChildValue(root, "Task", pTaskType, child=pTask, default='uploadcontacts')
    call getNodeName(pTaskType, buffer)

    call getChild(root, "Device", pDevice)
    call getChildValue(pDevice, "AtomRange", transpar%idxdevice)
    call getChild(pDevice, "FirstLayerAtoms", pTmp, requested=.false.)
    call readFirstLayerAtoms(pTmp, transpar%PL, transpar%nPLs, transpar%idxdevice)
    if (.not.associated(pTmp)) then
      call setChildValue(pDevice, "FirstLayerAtoms", transpar%PL)
    end if

    call getChildren(root, "Contact", pNodeList)
    transpar%ncont = getLength(pNodeList)
    allocate(transpar%contacts(transpar%ncont))
    call readContacts(pNodeList, transpar%contacts, geom, char(buffer), transpar%contactLayerTol)

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
        call detailedError(pTmp, strTmp)
      end if
      atomInRegion(transpar%contacts(ii)%idxrange(1):transpar%contacts(ii)%idxrange(2)) = .true.
    end do
    if (any(.not.atomInRegion)) then
      write(strTmp, "(A,I0,A)")"Atom ", findloc(atomInRegion, .false.),&
          & " is not in the device region or any contact"
      call detailedError(root, strTmp)
    end if

    call destroyNodeList(pNodeList)

    transpar%taskUpload = .false.

    select case (char(buffer))
    case ("contacthamiltonian")

      call getChildValue(pTaskType, "ContactId", buffer, child=pTmp)
      contact = getContactByName(transpar%contacts(:)%name, tolower(trim(unquote(char(buffer)))),&
          & pTmp)
      transpar%taskContInd = contact
      if (.not. geom%tPeriodic) then
        call getChildValue(pTaskType, "ContactSeparation", lateralContactSeparation, 1000.0_dp,&
            & modifier=modifier, child=field)
        call convertUnitHsd(char(modifier),lengthUnits,field,lateralContactSeparation)
      end if

      call reduceGeometry(transpar%contacts(contact)%lattice, transpar%contacts(contact)%idxrange,&
          & lateralContactSeparation, geom)

      transpar%ncont = 0

      call getChildValue(root, "writeBinaryContact", transpar%tWriteBinShift, .true.)

    case ("uploadcontacts")

      transpar%taskUpload = .true.

      call getChildValue(root, "readBinaryContact", transpar%tReadBinShift, .true.)

    case default

      call getNodeHSDName(pTaskType, buffer)
      call detailedError(pTask, "Invalid task '" // char(buffer) // "'")

   end select

   call destroyNodeList(pNodeList)

  end subroutine readTransportGeometry


  !> Reduce the geometry for the contact calculation
  subroutine reduceGeometry(contactVec, contactRange, lateralContactSeparation, geom)

    !> Vector between principle layers in the contact
    real(dp), intent(in) :: contactVec(3)

    !> Range of atoms in the contact
    integer, intent(in) :: contactRange(2)

    !> Lateral separation distance between contacts in a periodic box
    real(dp), intent(in) :: lateralContactSeparation

    !> atomic geometry
    type(TGeometry), intent(inout) :: geom

    real(dp) :: contUnitVec(3), dots(3), newLatVecs(3, 3), newOrigin(3)
    real(dp) :: minProj, maxProj
    logical :: mask(3)
    integer :: ind, indPrev, indNext, ii

    if (geom%tPeriodic) then
      contUnitVec = contactVec / sqrt(sum(contactVec**2, dim=1))
      dots = abs(matmul(contUnitVec, geom%latVecs))
      mask = (abs(dots - sqrt(sum(geom%latVecs, dim=1)**2)) < 1e-8_dp)
      if (count(mask) /= 1) then
        call error("Too many lattice vectors parallel to the contact")
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
  subroutine readFirstLayerAtoms(pnode, pls, npl, idxdevice, check)

    type(fnode), pointer, intent(in) :: pnode

    !> Start atoms in the principal layers
    integer, allocatable, intent(out) :: pls(:)

    !> Number of principal layers
    integer, intent(out) :: npl

    !> Atoms range of the device
    integer, intent(in) :: idxdevice(2)

    !> Optional setting to turn on/off check (defaults to on if absent)
    logical, optional, intent(in) :: check


    type(TListInt) :: li
    logical :: checkidx

    checkidx = .true.
    if (present(check)) checkidx = check

    if (associated(pnode)) then
        call init(li)
        call getChildValue(pnode, "", li)
        npl = len(li)
        allocate(pls(npl))
        call asArray(li, pls)
        call destruct(li)
        if (checkidx) then
          if (any(pls < idxdevice(1) .or. &
                  pls > idxdevice(2))) then
             call detailedError(pnode, "First layer atoms must be between " &
               &// i2c(idxdevice(1)) // " and " // i2c(idxdevice(2)) // ".")
          end if
        end if
      else
         npl = 1
         allocate(pls(npl))
         pls = (/ 1 /)
      end if

  end subroutine readFirstLayerAtoms


  !> Reads Green's function settings
  subroutine readGreensFunction(pNode, greendens, transpar, tempElec)

    !> Input tree
    type(fnode), pointer :: pTmp

    !> Settings for Green's function solver
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

    !> Transport solver settings
    type(TTransPar), intent(inout) :: transpar

    !> Electron temperature
    real(dp), intent(in) :: tempElec

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
      call getChildValue(pNode, "FermiLevel", fermiBuffer, modifier=modifier)
      if ( len(fermiBuffer) .eq. 1) then
        call asArray(fermiBuffer, greendens%oneFermi)
        greendens%oneFermi(2) = greendens%oneFermi(1)
      else if ( len(fermiBuffer) .eq. 2) then
        call asArray(fermiBuffer, greendens%oneFermi)
      else
        call detailedError(pNode, "FermiLevel accepts 1 or 2 (for collinear spin) values")
      end if
      call destruct(fermiBuffer)
      call convertUnitHsd(char(modifier), energyUnits, pNode, greendens%oneFermi)

      call getChild(pNode, "FirstLayerAtoms", pTmp, requested=.false.)
      call readFirstLayerAtoms(pTmp, greendens%PL, greendens%nPLs,&
                                &transpar%idxdevice, check = .false.)
      if (.not.associated(pTmp)) then
        call setChildValue(pNode, "FirstLayerAtoms", greendens%PL)
      end if
      !call getChild(pNode, "ContactPLs", pTmp, requested=.false.)
      !if (associated(pTmp)) then
      !  call init(li)
      !  call getChildValue(pTmp, "", li)
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

    call getChildValue(pNode, "LocalCurrents", greendens%doLocalCurr, .false.)
    call getChildValue(pNode, "Verbosity", greendens%verbose, 51)
    call getChildValue(pNode, "Delta", greendens%delta, 1.0e-5_dp, modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), energyUnits, field, greendens%delta)
    call getChildValue(pNode, "ReadSurfaceGFs", greendens%readSGF, .false.)
    call getChildValue(pNode, "SaveSurfaceGFs", greendens%saveSGF, .not.greendens%readSGF)
    call getChildValue(pNode, "ContourPoints", greendens%nP(1:2), [ 20, 20 ])
    call getChildValue(pNode, "EnclosedPoles",  greendens%nPoles, 3)
    call getChildValue(pNode, "LowestEnergy", greendens%enLow, -2.0_dp, modifier=modifier,&
        & child=field)
    call convertUnitHsd(char(modifier), energyUnits, field, greendens%enLow)
    call getChildValue(pNode, "FermiCutoff", greendens%nkT, 10)
      ! Fermi energy had not been set by other means yet

      ! Non equilibrium integration along real axis:
      ! The code will perform the integration if the number of points is larger
      ! than zero, no matter if there's bias or not.
      ! Therefore I restored the default on the energy step, as it works at zero
      ! bias and it scales flawlessy with increasing bias
      ! It is still allowed to directly set the number of points, if preferred
      ! libNEGF only wants the number of points in input
      call getChild(pNode, "RealAxisPoints", child1, requested=.false.)
      call getChild(pNode, "RealAxisStep", child2, requested=.false., &
          & modifier=buffer)
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
        call detailedError(child1, "RealAxisPoints and RealAxisStep " &
                            &// " cannot be specified together.")
      ! If only one is specified, take it as valid value
      else if (associated(child1)) then
        call getChildValue(pNode, "RealAxisPoints", greendens%nP(3))
      else if (associated(child2)) then
        call getChildValue(pNode, "RealAxisStep", Estep, child=child2, &
             & modifier=modifier)
        call convertUnitHsd(char(modifier), energyUnits, child2, Estep)
        realAxisConv = .true.
      ! If the system is under equilibrium we set the number of
      ! points to zero
      else if (equilibrium) then
        call getChildValue(pNode, "RealAxisPoints", greendens%nP(3), &
          & 0, child=child1)
      else
        !Default is a point every 1500H
        call getChildValue(pNode, "RealAxisStep", Estep, 6.65e-4_dp, &
                          &modifier=modifier, child=child2)
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
        !call getChildValue(pNode, "RealAxisPoints", greendens%nP(3), &
        !    & defvalue, child=child1)
      end if

  end subroutine readGreensFunction
#:endif


#:if WITH_POISSON

  !> Read in Poisson related data
#:if WITH_TRANSPORT
  subroutine readPoisson(pNode, poisson, tPeriodic, transpar, latVecs, updateSccAfterDiag)
#:else
  subroutine readPoisson(pNode, poisson, tPeriodic, latVecs, updateSccAfterDiag)
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
            call getChildValue(pNode, "PoissonThickness", poisson%poissBox(ii), modifier=modifier,&
                & child=field)
            call convertUnitHsd(char(modifier), lengthUnits, field, poisson%poissBox)
          else
            poisson%poissBox(ii) = sqrt(sum(latVecs(:,ii)**2))
          end if
        end do
      else
        call getChildValue(pNode, "PoissonBox", poisson%poissBox, modifier=modifier, child=field)
        call convertUnitHsd(char(modifier), lengthUnits, field, poisson%poissBox)
      end if
    #:else
      call getChildValue(pNode, "PoissonBox", poisson%poissBox, modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%poissBox)
    #:endif
    end if

    poisson%foundBox = needsPoissonBox
    call getChildValue(pNode, "MinimalGrid", poisson%poissGrid, [ 0.3_dp, 0.3_dp, 0.3_dp ],&
        & modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), lengthUnits, field, poisson%poissGrid)
    call getChildValue(pNode, "NumericalNorm", poisson%numericNorm, .false.)
    call getChild(pNode, "AtomDensityCutoff", pTmp, requested=.false., modifier=modifier)
    call getChild(pNode, "AtomDensityTolerance", pTmp2, requested=.false.)
    if (associated(pTmp) .and. associated(pTmp2)) then
      call detailedError(pNode, "Either one of the tags AtomDensityCutoff or AtomDensityTolerance&
          & can be specified.")
    else if (associated(pTmp)) then
      call getChildValue(pTmp, "", poisson%maxRadAtomDens, default=14.0_dp, modifier=modifier)
      call convertUnitHsd(char(modifier), lengthUnits, pTmp, poisson%maxRadAtomDens)
      if (poisson%maxRadAtomDens <= 0.0_dp) then
        call detailedError(pTmp2, "Atom density cutoff must be > 0")
      end if
    else
      call getChildValue(pNode, "AtomDensityTolerance", denstol, 1e-6_dp, child=pTmp2)
      if (denstol <= 0.0_dp) then
        call detailedError(pTmp2, "Atom density tolerance must be > 0")
      end if
      ! Negative value to signal automatic determination
      poisson%maxRadAtomDens = -denstol
    end if

    call getChildValue(pNode, "CutoffCheck", poisson%cutoffcheck, .true.)
    call getChildValue(pNode, "Verbosity", poisson%verbose, 51)
    call getChildValue(pNode, "SavePotential", poisson%savePotential, .false.)
    call getChildValue(pNode, "PoissonAccuracy", poisson%poissAcc, 1.0e-6_dp)
    call getChildValue(pNode, "BuildBulkPotential", poisson%bulkBC, .true.)
    call getChildValue(pNode, "ReadOldBulkPotential", poisson%readBulkPot, .false.)
    call getChildValue(pNode, "RecomputeAfterDensity", updateSccAfterDiag, .false.)
    call getChildValue(pNode, "MaxPoissonIterations", poisson%maxPoissIter, 60)

    poisson%overrideBC(:) = 0
    call getChild(pNode, "OverrideDefaultBC", pTmp, requested=.false.)
    if (associated(pTmp)) then
      call getPoissonBoundaryConditionOverrides(pTmp, [ 1, 2 ], poisson%overrideBC)
    end if

    call getChildValue(pNode, "OverrideBulkBC", pTmp, "none")
    poisson%overrBulkBC(:) = -1
    if (associated(pTmp)) then
      call getPoissonBoundaryConditionOverrides(pTmp, [ 0, 1, 2 ], poisson%overrBulkBC)
    end if

    call getChildValue(pNode, "BoundaryRegion", pTmp, "global")
    call getNodeName(pTmp, buffer)
    select case(char(buffer))
    case ("global")
      poisson%localBCType = "G"
    case ("square")
      poisson%localBCType = "S"
      call getChildValue(pTmp, "BufferLength", poisson%bufferLocBC, 9.0_dp, modifier=modifier,&
          & child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%bufferLocBC)
    case ("circle")
      poisson%localBCType = "C"
      call getChildValue(pTmp, "BufferLength", poisson%bufferLocBC, 9.0_dp, modifier=modifier,&
          & child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%bufferLocBC)
    case default
      call getNodeHSDName(pTmp, buffer)
      call detailedError(pTmp, "Invalid boundary region type '" // char(buffer) // "'")
    end select

    call getChildValue(pNode, "BoxExtension", poisson%bufferBox, 0.0_dp, modifier=modifier,&
        & child=field)
    call convertUnitHsd(char(modifier), lengthUnits, field, poisson%bufferBox)
    if (poisson%bufferBox.lt.0.0_dp) then
      call detailedError(pNode, "BoxExtension must be a positive number")
    endif

    ! PARSE GATE OPTIONS
    call getChildValue(pNode,"Gate",pTmp2,"none",child=pChild)
    call getNodeName(pTmp2, buffer)

    poisson%insLength = 0.0_dp
    poisson%insRad = 0.0_dp
    select case(char(buffer))
    case ("none")
      poisson%gateType = "N"
    case ("planar")
      poisson%gateType = "P"
      call getChildValue(pTmp2, "GateLength", poisson%gateLength_l, 0.0_dp, modifier= modifier,&
          & child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateLength_l)

      gatelength_l = poisson%gateLength_l !avoids a warning on intents
      call getChildValue(pTmp2, "GateLength_l", poisson%gateLength_l, gateLength_l,&
          & modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateLength_l)

      call getChildValue(pTmp2, "GateLength_t", poisson%gateLength_t, poisson%gateLength_l,&
          & modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateLength_t)

      call getChildValue(pTmp2, "GateDistance", poisson%gateRad, 0.0_dp, modifier=modifier,&
          & child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateRad)

      call getChildValue(pTmp2, "GatePotential", poisson%gatepot, 0.0_dp, modifier=modifier,&
          & child=field)
      call convertUnitHsd(char(modifier), energyUnits, field, poisson%gatepot)

      !call getChildValue(pTmp2, "GateDirection", poisson%gatedir, 2)
      poisson%gatedir = 2

    case ("cylindrical")
      poisson%gateType = "C"
      call getChildValue(pTmp2, "GateLength",poisson%gateLength_l, 0.0_dp, modifier= modifier,&
          & child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateLength_l)

      call getChildValue(pTmp2, "GateRadius", poisson%gateRad, 0.0_dp, modifier=modifier,&
          & child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, poisson%gateRad)

      call getChildValue(pTmp2, "GatePotential", poisson%gatepot, 0.0_dp, modifier=modifier,&
          & child=field)
      call convertUnitHsd(char(modifier), energyUnits, field, poisson%gatepot)

    case default
      call getNodeHSDName(pTmp2, buffer)
      call detailedError(pTmp2, "Invalid gate type '" // char(buffer) // "'")

    end select

    call getChildValue(pNode, "MaxParallelNodes", poisson%maxNumNodes, 1)

    poisson%scratch = "contacts"

  end subroutine readPoisson


  !> Over-rides the boundary conditions on the Poisson solver
  subroutine getPoissonBoundaryConditionOverrides(pNode, availableConditions, overrideBC)

    !> Input data tree
    type(fnode), pointer, intent(in) :: pNode

    !> List of conditions that can be set as choices
    integer, intent(in) :: availableConditions(:)

    !> Array of boundary condition types on the 6 faces of the box, 0 for use of default
    integer, intent(inout) :: overrideBC(:)

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
      call getChild(pNode, trim(bcstr(bctype)), pNode2, requested=.false.)
      if (associated(pNode2)) then
        call init(lStr)
        call getChildValue(pNode2, "boundaries", lStr, child=pChild)
        if (len(lStr).gt.6) then
          call detailedError(pChild,"boundaries must be 6 or less")
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
        call detailedError(pChild, &
            & "periodic override must be set both min max")
      end if
    end do

  end subroutine getPoissonBoundaryConditionOverrides

#:endif


#:if WITH_TRANSPORT
  !> Sanity checking of atom ranges and returning contact vector and direction.
  subroutine getContactVector(atomrange, geom, id, name, pContact, contactLayerTol, contactVec,&
      & contactDir)

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

    integer :: iStart, iStart2, iEnd, ii
    logical :: mask(3)
    character(lc) :: errorStr

    !! Sanity check for the atom ranges
    iStart = atomrange(1)
    iEnd = atomrange(2)
    if (iStart < 1 .or. iEnd < 1 .or. iStart > geom%nAtom .or. iEnd > geom%nAtom) then
      call detailedError(pContact, "Invalid atom range '" // i2c(iStart) &
          &// " " // i2c(iEnd) // "', values should be between " // i2c(1) &
          &// " and " // i2c(geom%nAtom) // ".")
    end if
    if (iEnd < iStart) then
      call detailedError(pContact, "Invalid atom order in contact '" // i2c(iStart) // " " //&
          & i2c(iEnd) // "', should be asscending order.")
    end if

    if (mod(iEnd - iStart + 1, 2) /= 0) then
      call detailedError(pContact, "Nr. of atoms in the contact must be even")
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
      call error(errorStr)

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
  subroutine readDephasing(node, orb, geom, tp, tundos)

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

    type(fnode), pointer :: value1, child

    call getChild(node, "VibronicElastic", child, requested=.false.)
    if (associated(child)) then
      tp%tDephasingVE = .true.
      call readElPh(child, tundos%elph, geom, orb, tp)
    end if

    call getChildValue(node, "BuettikerProbes", value1, "", child=child, &
        &allowEmptyValue=.true., dummyValue=.true.)
    if (associated(value1)) then
      tp%tDephasingBP = .true.
      call readDephasingBP(child, tundos%bp, geom, orb, tp)
    end if

    ! Lowdin transformations involve dense matrices and works only in small systems
    ! For the dftb+ official release the options are disabled
    tp%tOrthonormal = .false.
    tp%tOrthonormalDevice = .false.
    !call getChildValue(node, "Orthonormal", tp%tOrthonormal, .false.)
    !call getChildValue(node, "OrthonormalDevice", tp%tOrthonormalDevice, .false.)
    tp%tNoGeometry = .false.
    tp%NumStates = 0

  end subroutine readDephasing


  !> Read Electron-Phonon blocks (for density and/or current calculation)
  subroutine readElPh(node, elph, geom, orb, tp)

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


    logical :: block_model, semilocal_model

    elph%defined = .true.
    !! Only local el-ph model is defined (elastic for now)
    elph%model = 1

    call getChildValue(node, "MaxSCBAIterations", elph%scba_niter, default=100)
    call getChildValue(node, "atomBlock", block_model, default=.false.)
    if (block_model) then
      elph%model = 2
    endif

    !BUG: semilocal model crashes because of access of S before its allocation
    !     this because initDephasing was moved into initprogram
    call getChildValue(node, "semiLocal", semilocal_model, default=.false.)
    if (semilocal_model) then
      call detailedError(node, "semilocal dephasing causes crash and has been "//&
           & "temporarily disabled")
      elph%model = 3
    endif

    call readCoupling(node, elph, geom, orb, tp)

  end subroutine readElPh


  !> Read Buettiker probe dephasing blocks (for density and/or current calculation)
  subroutine readDephasingBP(node, elph, geom, orb, tp)

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

    logical :: block_model, semilocal_model
    type(string) :: model
    type(fnode), pointer :: dephModel

    call detailedError(node,"Buettiker probes are still under development")

    elph%defined = .true.
    call getChildValue(node, "", dephModel)
    call getNodeName2(dephModel, model)

    select case(char(model))
    case("dephasingprobes")
      !! Currently only zeroCurrent condition is implemented
      !! This corresponds to elastic dephasing probes
      tp%tZeroCurrent=.true.
      !! Only local bp model is defined (elastic for now)
    case("voltageprobes")
      call detailedError(dephModel,"voltageProbes have been not implemented yet")
      tp%tZeroCurrent=.false.
    case default
      call detailedError(dephModel,"unknown model")
    end select

    elph%model = 1

    call getChildValue(dephModel, "MaxSCBAIterations", elph%scba_niter, default=100)

    call getChildValue(dephModel, "atomBlock", block_model, default=.false.)
    if (block_model) then
      elph%model = 2
    endif

    !BUG: semilocal model crashes because of access of S before its allocation
    !     this because initDephasing occurs in initprogram
    call getChildValue(dephModel, "semiLocal", semilocal_model, default=.false.)
    if (semilocal_model) then
      call detailedError(dephModel, "semilocal dephasing is not working yet")
      elph%model = 3
    endif

    call readCoupling(dephModel, elph, geom, orb, tp)

  end subroutine readDephasingBP


  !> Reads coupling strength and mode for dephasing
  !> 2 modes support, constant or specified per each orbital
  subroutine readCoupling(node, elph, geom, orb, tp)

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

    call getChildValue(node, "Coupling", val, "", child=child, &
        & allowEmptyValue=.true., modifier=modifier, dummyValue=.true., list=.false.)

    call getNodeName(val, method)

    ! This reads also things like:  "Coupling [eV] = 0.34"
    !if (is_numeric(char(method))) then
    !  call getChildValue(node, "Coupling", rTmp, child=field)
    !  call convertUnitHsd(char(modifier), energyUnits, field, rTmp)
    !  elph%coupling = rTmp
    !  return
    !end if

    select case (char(method))
    case ("allorbitals")
      call getChild(child, "AllOrbitals", child2, requested=.false.)
      call getChildValue(child2, "", elph%coupling, child=field)
      call convertUnitHsd(char(modifier), energyUnits, field, elph%coupling)

    case ("atomcoupling")
      call getChild(child, "AtomCoupling", child2, requested=.false.)
      allocate(atmCoupling(atm_range(2)-atm_range(1)+1))
      atmCoupling = 0.d0
      call getChildren(child2, "AtomList", children)
      do ii = 1, getLength(children)
        call getItem1(children, ii, child3)
        call getChildValue(child3, "Atoms", buffer, child=child4, multiple=.true.)
        call getSelectedAtomIndices(child4, char(buffer), geom%speciesNames, geom%species, tmpI1)
        call getChildValue(child3, "Value", rTmp, child=field, modifier=modifier2)
        ! If not defined, use common unit modifier defined after Coupling
        if (len(modifier2)==0) then
          call convertUnitHsd(char(modifier), energyUnits, field, rTmp)
        else
          call convertUnitHsd(char(modifier2), energyUnits, field, rTmp)
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
      call getChildValue(child, "Constant", rtmp, child=field)
      call convertUnitHsd(char(modifier), energyUnits, field, rTmp)
      elph%coupling = rTmp

    case default
      call detailedError(node, "Coupling definition unknown")
    end select

  end subroutine readCoupling


  !> Read Tunneling and Dos options from analysis block
  subroutine readTunAndDos(root, orb, geo, tundos, transpar, tempElec)
    type(fnode), pointer :: root
    type(TOrbitals), intent(in) :: orb
    type(TGeometry), intent(in) :: geo

    !> tundos is the container to be filled
    type(TNEGFTunDos), intent(inout) :: tundos
    type(TTransPar), intent(inout) :: transpar
    real(dp), intent(in) :: tempElec

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

    call getChildValue(root, "Verbosity", tundos%verbose, 51)
    call getChildValue(root, "WriteLDOS", tundos%writeLDOS, .true.)
    call getChildValue(root, "WriteTunn", tundos%writeTunn, .true.)

    ! Read Temperature. Can override contact definition
    allocate(tundos%kbT(ncont))
    call getChild(root, "ContactTemperature", pTmp, modifier=modifier, requested=.false.)
    if (associated(pTmp)) then
      call init(temperature)
      call getChildValue(pTmp, "", temperature)
      if (len(temperature) .ne. ncont) then
        call detailedError(root, "ContactTemperature does not match the number of contacts")
      end if
      call asArray(temperature, tundos%kbT)
      call destruct(temperature)
      call convertUnitHsd(char(modifier), energyUnits, pTmp, tundos%kbT)
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
      call getChildValue(root, "EnergyRange", eRange, modifier=modifier,&
      & child=field)
      call convertUnitHsd(char(modifier), energyUnits, field, eRange)
      call getChildValue(root, "EnergyStep", tundos%estep,&
      & modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), energyUnits, field, tundos%estep)
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
      call getChildValue(root, "EnergyStep", tundos%estep, 6.65e-4_dp, &
                          &modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), energyUnits, field, tundos%estep)
      call getChildValue(root, "EnergyRange", eRange, eRangeDefault, &
                          modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), energyUnits, field, eRange)
    end if

    tundos%emin = eRange(1)
    tundos%emax = eRange(2)
    ! Terminal currents
    call getChild(root, "TerminalCurrents", pTmp, requested=.false.)
      if (associated(pTmp)) then
        call getChildren(pTmp, "EmitterCollector", pNodeList)
        allocate(tundos%ni(getLength(pNodeList)))
        allocate(tundos%nf(getLength(pNodeList)))
        do ii = 1, getLength(pNodeList)
          call getItem1(pNodeList, ii, pNode)
          call getEmitterCollectorByName(pNode, tundos%ni(ii), tundos%nf(ii),&
              & transpar%contacts(:)%name)
        end do
        call destroyNodeList(pNodeList)
      else
        allocate(tundos%ni(ncont-1) )
        allocate(tundos%nf(ncont-1) )
        call setChild(root, "TerminalCurrents", pTmp)
        ind = 1
        do ii = 1, 1
          do jj = ii + 1, ncont
            call setChildValue(pTmp, "EmitterCollector", &
                &(/ transpar%contacts(ii)%name, transpar%contacts(jj)%name /))
            tundos%ni(ind) = ii
            tundos%nf(ind) = jj
            ind = ind + 1
          end do
        end do
      end if
      call getChildValue(root, "Delta", tundos%delta, &
          &1.0e-5_dp, modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), energyUnits, field, &
          &tundos%delta)
      call getChildValue(root, "BroadeningDelta", tundos%broadeningDelta, &
          &0.0_dp, modifier=modifier, child=field)
      call convertUnitHsd(char(modifier), energyUnits, field, &
          &tundos%broadeningDelta)

      call readPDOSRegions(root, geo, transpar%idxdevice, iAtInRegion, &
          & tShellResInRegion, regionLabelPrefixes)

      if (allocated(iAtInRegion)) then
        call transformPdosRegionInfo(iAtInRegion, tShellResInRegion, &
            & regionLabelPrefixes, orb, geo%species, tundos%dosOrbitals, &
            & tundos%dosLabels)
      end if

  end subroutine readTunAndDos


  !> Read bias information, used in Analysis and Green's function eigensolver
  subroutine readContacts(pNodeList, contacts, geom, task, contactLayerTol)

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

    integer :: ii
    type(fnode), pointer :: field, pNode, pTmp, child1, child2
    type(string) :: buffer, modifier
    type(TListReal) :: fermiBuffer

    do ii = 1, size(contacts)

      contacts(ii)%wideBand = .false.
      contacts(ii)%wideBandDos = 0.0_dp

      call getItem1(pNodeList, ii, pNode)
      call getChildValue(pNode, "Id", buffer, child=pTmp)
      buffer = tolower(trim(unquote(char(buffer))))
      if (len(buffer) > mc) then
        call detailedError(pTmp, "Contact id may not be longer than " // i2c(mc) // " characters.")
      end if
      contacts(ii)%name = char(buffer)
      if (any(contacts(1:ii-1)%name == contacts(ii)%name)) then
        call detailedError(pTmp, "Contact id '" // trim(contacts(ii)%name) //  "' already in use")
      end if

      call getChildValue(pNode, "PLShiftTolerance", contactLayerTol, 1e-5_dp, modifier=modifier,&
          & child=field)
      call convertUnitHsd(char(modifier), lengthUnits, field, contactLayerTol)

      call getChildValue(pNode, "AtomRange", contacts(ii)%idxrange, child=pTmp)
      call getContactVector(contacts(ii)%idxrange, geom, ii, contacts(ii)%name, pTmp,&
        & contactLayerTol, contacts(ii)%lattice, contacts(ii)%dir)
      contacts(ii)%length = sqrt(sum(contacts(ii)%lattice**2))

      ! Contact temperatures. A negative default is used so it is quite clear when the user sets a
      ! different value. In such a case this overrides values defined in the Filling block
      call getChild(pNode,"Temperature", field, modifier=modifier, requested=.false.)
      if (associated(field)) then
        call getChildValue(pNode, "Temperature", contacts(ii)%kbT, 0.0_dp, modifier=modifier,&
            & child=field)
        call convertUnitHsd(char(modifier), energyUnits, field, contacts(ii)%kbT)
      else
        contacts(ii)%kbT = -1.0_dp ! -1.0 simply means 'not defined'
      end if

      if (task .eq. "uploadcontacts") then
        call getChildValue(pNode, "Potential", contacts(ii)%potential, 0.0_dp, modifier=modifier,&
            & child=field)
        call convertUnitHsd(char(modifier), energyUnits, field, contacts(ii)%potential)

        call getChildValue(pNode, "WideBand", contacts(ii)%wideBand, .false.)

        if (contacts(ii)%wideBand) then

          ! WideBandApproximation is defined as energy spacing between levels of the contact. In the
          ! code the inverse value (Density of states) is used. Convert the negf input
          ! value. Default is 20 / e eV.
          call getChildValue(pNode, "LevelSpacing", contacts(ii)%wideBandDos, 0.735_dp,&
              & modifier=modifier, child=field)
          call convertUnitHsd(char(modifier), energyUnits, field, contacts(ii)%wideBandDos)
          contacts(ii)%wideBandDos = 1.d0 / contacts(ii)%wideBandDos

        end if


        ! Fermi level: in case of collinear spin we accept two values (up and down)
        ! call init(fermiBuffer)
        ! call getChildValue(pNode, "FermiLevel", fermiBuffer, modifier=modifier)
        ! if ( len(fermiBuffer) .eq. 1) then
        !   call asArray(fermiBuffer, contacts(ii)%eFermi)
        !   contacts(ii)%eFermi(2) = contacts(ii)%eFermi(1)
        ! else if ( len(fermiBuffer) .eq. 2) then
        !   call asArray(fermiBuffer, contacts(ii)%eFermi)
        ! else
        !   call detailedError(pNode, "FermiLevel accepts 1 or 2 (for collinear spin) values")
        ! end if
        ! call destruct(fermiBuffer)


        call getChildValue(pNode, "FermiLevel", child1, "", child=child2, allowEmptyValue=.true.,&
            & modifier=modifier)
        call getNodeName2(child1, buffer)
        if (char(buffer) == "") then
          contacts(ii)%tFermiSet = .false.
          call detailedWarning(pNode, "Missing Fermi level - required to be set in solver block or&
              & read from a contact shift file")
        else
          call init(fermiBuffer)
          call getChildValue(child2, "", fermiBuffer, modifier=modifier)
          select case(len(fermiBuffer))
          case (1)
            call asArray(fermiBuffer, contacts(ii)%eFermi(1:1))
            contacts(ii)%eFermi(2) = contacts(ii)%eFermi(1)
          case (2)
            call asArray(fermiBuffer, contacts(ii)%eFermi(:2))
          case default
            call detailedError(pNode, "FermiLevel accepts 1 or 2 (for collinear spin) values")
          end select
          call destruct(fermiBuffer)
          call convertUnitHsd(char(modifier), energyUnits, child2, contacts(ii)%eFermi)

          contacts(ii)%tFermiSet = .true.

          ! NOTE: These options have been commented out: there is a problem in parallel execution
          ! since one single file is accessed by all processors causing rush conditions
          ! The options are therefore disabled for the official dftb+ release
          !call getChildValue(pNode, "WriteSelfEnergy", contacts(ii)%tWriteSelfEnergy, .false.)
          !call getChildValue(pNode, "WriteSurfaceGF", contacts(ii)%tWriteSurfaceGF, .false.)
          !call getChildValue(pNode, "ReadSelfEnergy", contacts(ii)%tReadSelfEnergy, .false.)
          !call getChildValue(pNode, "ReadSurfaceGF", contacts(ii)%tReadSurfaceGF, .false.)
          contacts(ii)%tWriteSelfEnergy = .false.
          contacts(ii)%tWriteSurfaceGF = .false.
          contacts(ii)%tReadSelfEnergy = .false.
          contacts(ii)%tReadSurfaceGF = .false.
        end if

      end if

    end do

  end subroutine readContacts


  !> Read in Fermi levels
  subroutine getFermiLevels(pNode, eFermis, nodeModifier)

    !> Document tree node to start from
    type(fnode), pointer :: pNode

    !> Fermi energies for contacts
    real(dp), intent(out) :: eFermis(:)

    !> Any node modifiers in action
    type(string), intent(in) :: nodeModifier

    real(dp) :: eFermi
    type(fnode), pointer :: pChild
    type(string) :: modifier

    call getChild(pNode, "SetForAll", pChild, requested=.false.)
    if (associated(pChild)) then
      call getChildValue(pChild, "", eFermi)
      call convertUnitHsd(char(nodeModifier), energyUnits, pNode, eFermi)
      eFermis(:) = eFermi
    else
      call getChildValue(pNode, "", eFermis, modifier=modifier, child=pChild)
      call convertUnitHsd(char(modifier), energyUnits, pChild, eFermis)
    end if

  end subroutine getFermiLevels


  !> Get contacts for terminal currents by name
  subroutine getEmitterCollectorByName(pNode, emitter, collector, contactNames)

    !> Node in the input tree for error reporting
    type(fnode), pointer :: pNode

    !> Contact number for emitting
    integer, intent(out) :: emitter

    !> Contact number for collecting
    integer, intent(out) :: collector

    !> Labels of contacts
    character(len=*), intent(in) :: contactNames(:)

    type(TListString) :: lString
    character(len=mc) :: buffer

    call init(lString)
    call getChildValue(pNode, "", lString)
    if (len(lString) /= 2) then
      call detailedError(pNode, "You must provide two contacts")
    end if
    call get(lString, buffer, 1)
    emitter = getContactByName(contactNames, buffer, pNode)
    call get(lString, buffer, 2)
    collector = getContactByName(contactNames, buffer, pNode)
    call destruct(lString)

  end subroutine getEmitterCollectorByName


  !> Getting the contact by name
  function getContactByName(contactNames, contName, pNode) result(contact)

    !> Node in the input tree for error reporting
    type(fnode), pointer :: pNode

    !> All of the contact labels
    character(len=*), intent(in) :: contactNames(:)

    !> Specific contact label to identify
    character(len=*), intent(in) :: contName

    !> Contact number
    integer :: contact

    logical :: tFound

    tFound = .false.
    do contact = 1, size(contactNames)
      tFound = (contactNames(contact) == contName)
      if (tFound) then
        exit
      end if
    end do
    if (.not. tFound) then
      call detailedError(pNode, "Invalid collector contact name '" // trim(contName) // "'")
    end if

  end function getContactByName


  !> Read the names of regions to calculate PDOS for
  subroutine readPDOSRegions(node, geo, idxdevice, iAtInregion, tShellResInRegion, regionLabels)

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
      call getChildValue(node, "ComputeLDOS", do_ldos, .true.)
      if (do_ldos) then
        write(strTmp,"(I0, ':', I0)") idxdevice(1), idxdevice(2)
        call setChild(node, "Region", child)
        call setChildValue(child, "Atoms", trim(strTmp))
        call setChildValue(child, "Label", "localDOS")
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
      call getChildValue(child, "Atoms", buffer, child=child2, multiple=.true.)
      call getSelectedAtomIndices(child2, char(buffer), geo%speciesNames,&
          & geo%species(idxdevice(1) : idxdevice(2)), tmpI1,&
          & selectionRange=[idxdevice(1), idxdevice(2)], indexRange=[1, geo%nAtom])
      iAtInRegion(iReg)%data = tmpI1
      call getChildValue(child, "ShellResolved", tShellResInRegion(iReg), .false., child=child2)
      if (tShellResInRegion(iReg)) then
        if (.not. all(geo%species(tmpI1) == geo%species(tmpI1(1)))) then
          call detailedError(child2, "Shell resolved PDOS can only summed up over atoms of the same&
              & type")
        end if
      end if
      write(strTmp, "('region',I0)") iReg
      call getChildValue(child, "Label", buffer, trim(strTmp))
      regionLabels(iReg) = unquote(char(buffer))
    end do
    call destroyNodeList(children)

  end subroutine readPDOSRegions


  !> Some assignment and consistency check in negf/poisson containers before calling initialization
  subroutine finalizeNegf(input)

    !> Input structure for DFTB+
    type(TInputData), intent(inout) :: input

    integer :: ii

    !! Check consistency between different deltas
    if (input%ginfo%tundos%defined.and.input%ginfo%greendens%defined) then
      if (input%ginfo%tundos%delta.ne.input%ginfo%greendens%delta) then
        call error("Delta parameter must be the same in GreensFunction and TunnelingAndDos")
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
            call error("Contact " // i2c(ii) // " not parallel to any &
              & coordinate axis is not compatible with Poisson solver")
          end if
        end do
      end if
    #:endblock

    !! Temporarily not supporting surface green function read/load
    !! for spin polarised, because spin is handled outside of libnegf
    if (input%ginfo%greendens%defined) then
      if (input%ctrl%tSpin .and. input%ginfo%greendens%saveSGF) then
        call error("SaveSurfaceGFs must be disabled in collinear spin calculations")
      end if
      if  (input%ctrl%tSpin .and. input%ginfo%greendens%readSGF) then
        call error("ReadSurfaceGFs must be disabled in collinear spin calculations")
      end if
    end if

  end subroutine finalizeNegf
#:endif


  !> This subroutine overrides the neutral (reference) atom electronic occupation
  subroutine readCustomReferenceOcc(root, orb, referenceOcc, geo, iAtInRegion, customOcc)

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

    type(fnode), pointer :: node, container, child
    type(fnodeList), pointer :: nodes
    type(string) :: buffer
    integer :: nCustomOcc, iCustomOcc, iShell, iSpecies, nAtom
    character(sc), allocatable :: shellNamesTmp(:)
    logical, allocatable :: atomOverriden(:)

    call renameChildren(root, "CustomizedOccupations", "CustomisedOccupations")
    call getChild(root, "CustomisedOccupations", container, requested=.false.)
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
      call getChildValue(node, "Atoms", buffer, child=child, multiple=.true.)
      call getSelectedAtomIndices(child, char(buffer), geo%speciesNames, geo%species,&
          & iAtInRegion(iCustomOcc)%data)
      if (any(atomOverriden(iAtInRegion(iCustomOcc)%data))) then
        call detailedError(child, "Atom region contains atom(s) which have already been overridden")
      end if
      atomOverriden(iAtInRegion(iCustomOcc)%data) = .true.
      iSpecies = geo%species(iAtInRegion(iCustomOcc)%data(1))
      if (any(geo%species(iAtInRegion(iCustomOcc)%data) /= iSpecies)) then
        call detailedError(child, "All atoms in a ReferenceOccupation declaration must have the&
            & same type.")
      end if
      call getShellNames(iSpecies, orb, shellNamesTmp)
      do iShell = 1, orb%nShell(iSpecies)
          call getChildValue(node, shellNamesTmp(iShell), customOcc(iShell, iCustomOcc), &
            & default=referenceOcc(iShell, iSpecies))
      end do
      deallocate(shellNamesTmp)
    end do
    call destroyNodeList(nodes)

  end subroutine readCustomReferenceOcc


  !> Reads the parallel block.
  subroutine readParallel(root, input)

    !> Root node eventually containing the current block
    type(fnode), pointer, intent(in) :: root

    !> Input structure to be filled
    type(TInputData), intent(inout) :: input

    type(fnode), pointer :: node, pTmp

    call getChild(root, "Parallel", child=node, requested=.false.)
    if (withMpi .and. .not. associated(node)) then
      call setChild(root, "Parallel", node)
    end if
    if (associated(node)) then
      if (.not. withMpi) then
        call detailedWarning(node, "Settings will be read but ignored (compiled without MPI&
            & support)")
      end if
      allocate(input%ctrl%parallelOpts)
      call getChildValue(node, "Groups", input%ctrl%parallelOpts%nGroup, 1, child=pTmp)
      call getChildValue(node, "UseOmpThreads", input%ctrl%parallelOpts%tOmpThreads, .not. withMpi)
      call readBlacs(node, input%ctrl%parallelOpts%blacsOpts)
    end if

  end subroutine readParallel


  !> Reads the blacs block
  subroutine readBlacs(root, blacsOpts)

    !> Root node eventually containing the current block
    type(fnode), pointer, intent(in) :: root

    !> Blacs settings
    type(TBlacsOpts), intent(inout) :: blacsOpts

    type(fnode), pointer :: node

    call getChild(root, "Blacs", child=node, requested=.false.)
    if (withScalapack .and. .not. associated(node)) then
      call setChild(root, "Blacs", node)
    end if
    if (associated(node)) then
      if (.not. withScalapack) then
        call detailedWarning(node, "Settings will be read but ignored (compiled without SCALAPACK&
            & support)")
      end if
      call getChildValue(node, "BlockSize", blacsOpts%blockSize, 32)
    end if

  end subroutine readBlacs


  !> Reads the settings for electrostatic potential plotting
  subroutine readElectrostaticPotential(node, geo, ctrl)

    !> Node containing optional electrostatic settings
    type(fnode), pointer, intent(in) :: node

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Control structure
    type(TControl), intent(inout) :: ctrl

    type(fnode), pointer :: child, child2, child3
    type(string) :: buffer, modifier
    type(TListRealR1) :: lr1

    call getChild(node, "ElectrostaticPotential", child, requested=.false.)
    if (.not. associated(child)) then
      return
    end if

    if (.not. ctrl%tSCC) then
      call error("Electrostatic potentials only available in an SCC calculation")
    end if
    allocate(ctrl%elStatPotentialsInp)
    call getChildValue(child, "OutputFile", buffer, "ESP.dat")
    ctrl%elStatPotentialsInp%espOutFile = unquote(char(buffer))
    ctrl%elStatPotentialsInp%tAppendEsp = .false.
    if (ctrl%isGeoOpt .or. ctrl%tMD) then
      call getChildValue(child, "AppendFile", ctrl%elStatPotentialsInp%tAppendEsp, .false.)
    end if
    call init(lr1)
    ! discrete points
    call getChildValue(child, "Points", child2, "", child=child3, modifier=modifier,&
        & allowEmptyValue=.true.)
    call getNodeName2(child2, buffer)
    if (char(buffer) /= "") then
      call getChildValue(child3, "", 3, lr1, modifier=modifier)
      allocate(ctrl%elStatPotentialsInp%espGrid(3,len(lr1)))
      call asArray(lr1, ctrl%elStatPotentialsInp%espGrid)
      if (geo%tPeriodic .and. (char(modifier) == "F" .or. char(modifier) == "f")) then
        ctrl%elStatPotentialsInp%espGrid = matmul(geo%latVecs, ctrl%elStatPotentialsInp%espGrid)
      else
        call convertUnitHsd(char(modifier), lengthUnits, child3,&
            & ctrl%elStatPotentialsInp%espGrid)
      end if
    end if
    call destruct(lr1)

    ! grid specification for points instead
    call getChild(child, "Grid", child=child2, modifier=modifier, requested=.false.)
    if (associated(child2)) then
      if (allocated(ctrl%elStatPotentialsInp%espGrid)) then
        call error("Both grid and point specification not both currently possible")
      end if
      if (geo%tPeriodic) then
        call readGrid(ctrl%elStatPotentialsInp%espGrid, child2, modifier,&
            & latVecs=geo%latVecs, nPoints=ctrl%elStatPotentialsInp%gridDimensioning,&
            & origin=ctrl%elStatPotentialsInp%origin,&
            & axes=ctrl%elStatPotentialsInp%axes)
      else
        call readGrid(ctrl%elStatPotentialsInp%espGrid, child2, modifier,&
            & nPoints=ctrl%elStatPotentialsInp%gridDimensioning,&
            & origin=ctrl%elStatPotentialsInp%origin,&
            & axes=ctrl%elStatPotentialsInp%axes)
      end if
    end if
    if (.not.allocated(ctrl%elStatPotentialsInp%espGrid)) then
      call detailedError(child,"Either a grid or set of points must be specified")
    end if
    call getChildValue(child, "Softening", ctrl%elStatPotentialsInp%softenESP, 1.0E-6_dp,&
        & modifier=modifier, child=child2)
    call convertUnitHsd(char(modifier), lengthUnits, child2, ctrl%elStatPotentialsInp%softenEsp)

  end subroutine readElectrostaticPotential


  !> Read in a grid specification
  subroutine readGrid(points, node, modifier, latVecs, nPoints, origin, axes)

    !> Points in the grid
    real(dp), allocatable, intent(out) :: points(:,:)

    !> input data to parse
    type(fnode), pointer, intent(in) :: node

    !> unit modifier for the grid
    type(string), intent(in) :: modifier

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
          & geometries")
    end if

    call getChildValue(node, "Spacing", r3Tmp, child=child)
    call getChildValue(node, "Origin", r3Tmpb, child=child)
    call getChildValue(node, "GridPoints", i3Tmp, child=child)
    if (any(i3Tmp < 1)) then
      call detailedError(child,"Grid must be at least 1x1x1")
    end if
    if (any(abs(r3Tmp) < epsilon(1.0_dp) .and. i3Tmp > 1)) then
      call detailedError(child,"Grid spacings must be non-zero")
    end if
    allocate(points(3,product(i3Tmp)))
    if (present(nPoints)) then
      nPoints = i3Tmp
    end if

    !  length not fraction modifier
    if (.not.(tPeriodic .and. (char(modifier) == "F" .or. char(modifier) == "f"))) then
      call convertUnitHsd(char(modifier), lengthUnits, child, r3Tmp)
      call convertUnitHsd(char(modifier), lengthUnits, child, r3Tmpb)
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
      call getChildValue(node, "Directions", axes_, r33Tmp, child=child)
      if (abs(determinant33(axes_)) < epsilon(1.0_dp)) then
        call detailedError(child, "Dependent axis directions")
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


  !> Parses hybrid xc-functional input.
  subroutine parseHybridBlock(node, input, geo, skFiles)

    !> Node to parse
    type(fnode), intent(in), pointer :: node

    !> Range separated data structure to fill
    type(THybridXcInp), intent(inout), allocatable :: input

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    !> List of SK file names to read in for every interaction
    type(TListCharLc), intent(inout) :: skFiles(:,:)

    !! File name of representative SK-file to read
    character(lc) :: fileName

    !! True, if hybrid xc-functional input block present
    logical :: isHybridInp

    !! True, if hybrid xc-functional extra tag found in SK-file(s)
    logical :: isHybridSk

    !! Hybrid xc-functional extra tag in SK-files, if allocated
    integer :: hybridXcSkTag

    !! Hybrid functional type of user input
    integer :: hybridXcInputTag

    !! Auxiliary node pointers
    type(fnode), pointer :: hybridChild, hybridValue, screeningChild
    type(fnode), pointer :: screeningValue, cmChild, cmValue, child1, child2

    !! Temporary string buffers
    type(string) :: buffer, modifier

    !! Temporary string buffer, that stores the gamma function type
    type(string) :: strBuffer
    character(lc) :: strTmp

    @:ASSERT(size(skFiles, dim=1) == size(skFiles, dim=2))
    @:ASSERT((size(skFiles, dim=1) > 0))

    ! Extracting hybridXc tag from first SK-file only is a workaround and assumes that a set of
    ! given SK-files uses the same parameters (which should always be the case)!
    call get(skFiles(1, 1), fileName, 1)

    ! Check if SK-files contain extra tag for hybrid xc-functionals
    call inquireHybridXcTag(fileName, hybridXcSkTag)
    isHybridSk = hybridXcSkTag /= hybridXcFunc%none

    call getChild(node, "Hybrid", child=hybridChild, requested=.false.)
    call getChildValue(node, "Hybrid", hybridValue, "None", child=hybridChild)
    call getNodeName(hybridValue, buffer)

    isHybridInp = associated(hybridChild) .and. (tolower(char(buffer)) /= "none")

    if (isHybridInp .and. .not. isHybridSk) then
      call error("Hybrid input block present, but SK-file '" // trim(fileName)&
          & // "' appears to be (semi-)local.")
    elseif (isHybridSk .and. .not. isHybridInp) then
      call error("Hybrid SK-file '" // trim(fileName) // "' present, but HSD input block missing.")
    end if

    if (isHybridInp) then
      ! Convert hybrid functional type of user input to enumerator
      select case(tolower(char(buffer)))
      case ("global")
        hybridXcInputTag = hybridXcFunc%hyb
      case ("lc")
        hybridXcInputTag = hybridXcFunc%lc
      case ("cam")
        hybridXcInputTag = hybridXcFunc%cam
      case default
        call detailedError(hybridChild, "Unknown hybrid xc-functional type '" // char(buffer)&
            & // "' in input.")
      end select

      ! Check if hybrid functional type is in line with SK-files
      if (.not. hybridXcInputTag == hybridXcSkTag) then
        call detailedError(hybridChild, "Hybrid functional type conflict with SK-files.")
      end if

      allocate(input)
      input%hybridXcType = hybridXcInputTag
      call getChildValue(hybridValue, "Screening", screeningValue, "Thresholded",&
          & child=screeningChild)

      call getNodeName(screeningValue, buffer)
      select case(tolower(char(buffer)))
      case ("neighbourbased")
        input%hybridXcAlg = hybridXcAlgo%neighbourBased
        call getChildValue(screeningValue, "CutoffReduction", input%cutoffRed, 0.0_dp,&
            & modifier=modifier, child=child1)
        call convertUnitHsd(char(modifier), lengthUnits, child1, input%cutoffRed)
        if (geo%tPeriodic) then
          call getChildValue(screeningValue, "Threshold", input%screeningThreshold, 1e-6_dp)
        end if
      case ("thresholded")
        input%hybridXcAlg = hybridXcAlgo%thresholdBased
        call getChildValue(screeningValue, "Threshold", input%screeningThreshold, 1e-6_dp)
        call getChildValue(screeningValue, "CutoffReduction", input%cutoffRed, 0.0_dp,&
            & modifier=modifier, child=child1)
        call convertUnitHsd(char(modifier), lengthUnits, child1, input%cutoffRed)
      case ("matrixbased")
        input%hybridXcAlg = hybridXcAlgo%matrixBased
        ! In this case, CutoffRedunction is not used so it should be set to zero.
        input%cutoffRed = 0.0_dp
      case default
        call getNodeHSdName(screeningValue, buffer)
        call detailedError(screeningChild, "Invalid screening method '" // char(buffer) // "'")
      end select

      ! Additional settings for periodic sytems
      ifPeriodic: if (geo%tPeriodic) then

        ! parse gamma function type (full, truncated, mic, ...)
        call getChildValue(hybridValue, "CoulombMatrix", cmValue, "Truncated", child=cmChild)

        call getNodeName(cmValue, buffer)
        select case(char(buffer))
        case ("full")
          input%gammaType = hybridXcGammaTypes%full
        case ("minimumimage")
          input%gammaType = hybridXcGammaTypes%mic
        case ("truncated")
          input%gammaType = hybridXcGammaTypes%truncated
        case ("truncated+damping")
          input%gammaType = hybridXcGammaTypes%truncatedAndDamped
        case default
          call getNodeHSdName(cmValue, buffer)
          call detailedError(cmChild, "Invalid Gamma function type '" // char(strBuffer) // "'")
        end select

        ! g-Summation cutoff not needed for MIC CAM Hamiltonian
        if (input%gammaType /= hybridXcGammaTypes%mic) then
          call getChild(cmValue, "GSummationCutoff", child=child1, modifier=modifier,&
              & requested=.false.)
          if (associated(child1)) then
            allocate(input%gSummationCutoff)
            call getChildValue(child1, "", input%gSummationCutoff, modifier=modifier, child=child2)
            call convertUnitHsd(char(modifier), lengthUnits, child2, input%gSummationCutoff)
          end if
        end if

        if (input%gammaType == hybridXcGammaTypes%truncated&
            & .or. input%gammaType == hybridXcGammaTypes%truncatedAndDamped) then
          call getChild(cmValue, "CoulombCutoff", child=child1, modifier=modifier,&
              & requested=.false.)
          if (associated(child1)) then
            allocate(input%gammaCutoff)
            call getChildValue(child1, "", input%gammaCutoff, modifier=modifier, child=child2)
            call convertUnitHsd(char(modifier), lengthUnits, child2, input%gammaCutoff)
          end if
        end if

      else
        ! Always use unaltered gamma function for non-periodic systems
        input%gammaType = hybridXcGammaTypes%full
      end if ifPeriodic

      ! Number of primitive cells regarded in MIC, along each supercell folding direction
      if (input%gammaType == hybridXcGammaTypes%mic) then
        allocate(input%wignerSeitzReduction)
        call getChildValue(cmValue, "WignerSeitzReduction", input%wignerSeitzReduction, default=0)
      end if

    else

      hybridXcInputTag = hybridXcFunc%none

    end if

  end subroutine parseHybridBlock


  !> Reads the REKS block
  subroutine readReks(node, dummy, ctrl, geo)

    !> Node to parse
    type(fnode), pointer, intent(in) :: node

    !> Node to parse
    type(fnode), pointer, intent(in) :: dummy

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    type(string) :: buffer

    ! SSR(2,2) or SSR(4,4) stuff
    call getNodeName(dummy, buffer)

    select case (char(buffer))
    case ("none")
      ctrl%reksInp%reksAlg = reksTypes%noReks
    case ("ssr22")
      ctrl%reksInp%reksAlg = reksTypes%ssr22
      call readSSR22(dummy, ctrl, geo)
    case ("ssr44")
      ctrl%reksInp%reksAlg = reksTypes%ssr44
      call detailedError(node, "SSR(4,4) is not implemented yet.")
    case default
      call getNodeHSDName(dummy, buffer)
      call detailedError(node, "Invalid Algorithm '" // char(buffer) // "'")
    end select

  end subroutine readReks


  !> Reads the SSR(2,2) block
  subroutine readSSR22(node, ctrl, geo)

    !> Node to parse
    type(fnode), pointer, intent(in) :: node

    !> Control structure to fill
    type(TControl), intent(inout) :: ctrl

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    type(fnode), pointer :: child1, value2, child2
    type(TListString) :: strBuffer
    type(string) :: buffer2
    character(sc), allocatable :: tmpFunc(:)
    integer :: ii, nFunc
    logical :: tFunc = .true.


    !> Read 'Energy' block
    call getChild(node, "Energy", child=child1)

    !> Read 'Functional' block in 'Energy' block
    call init(strBuffer)
    call getChildValue(child1, "Functional", strBuffer)
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
      call detailedError(child1, "Invalid Functional")
    end if

    !> Decide the energy states in SA-REKS
    !> If true, it includes all possible states in current active space
    !> If false, it includes the states used in minimized energy functional
    call getChildValue(child1, "IncludeAllStates", ctrl%reksInp%tAllStates, default=.false.)
    !> Calculate SSR state with inclusion of SI, otherwise calculate SA-REKS state
    call getChildValue(child1, "StateInteractions", ctrl%reksInp%tSSR, default=.false.)


    !> Target SSR state
    call getChildValue(node, "TargetState", ctrl%reksInp%rstate, default=1)
    !> Target microstate
    call getChildValue(node, "TargetMicrostate", ctrl%reksInp%Lstate, default=0)

    !> Read initial guess for eigenvectors in REKS
    !> If true, initial eigenvectors are obtained from 'eigenvec.bin'
    !> If false, initial eigenvectors are obtained from diagonalisation of H0
    call getChildValue(node, "ReadEigenvectors", ctrl%reksInp%tReadMO, default=.false.)
    !> Maximum iteration used in FON optimisation
    call getChildValue(node, "FonMaxIter", ctrl%reksInp%FonMaxIter, default=20)
    !> Shift value in SCC cycle
    call getChildValue(node, "Shift", ctrl%reksInp%shift, default=0.3_dp)

    !> Read "SpinTuning" block with 'nType' elements
    call readSpinTuning(node, ctrl, geo%nSpecies)

    !> Calculate transition dipole moments
    call getChildValue(node, "TransitionDipole", ctrl%reksInp%tTDP, default=.false.)


    !> Read 'Gradient' block
    !> Algorithms to calculate analytical gradients
    call getChildValue(node, "Gradient", value2, "ConjugateGradient", child=child2)
    call getNodeName(value2, buffer2)

    select case (char(buffer2))
    case ("conjugategradient")
      !> Maximum iteration used in calculation of gradient with PCG and CG
      call getChildValue(value2, "CGmaxIter", ctrl%reksInp%CGmaxIter, default=20)
      !> Tolerance used in calculation of gradient with PCG and CG
      call getChildValue(value2, "Tolerance", ctrl%reksInp%Glimit, default=1.0E-8_dp)
      !> Use preconditioner for conjugate gradient algorithm
      call getChildValue(value2, "Preconditioner", ctrl%reksInp%tPrecond, default=.false.)
      !> Save 'A' and 'Hxc' to memory in gradient calculation
      call getChildValue(value2, "SaveMemory", ctrl%reksInp%tSaveMem, default=.false.)
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
      call detailedError(child2, "Invalid Algorithm '" // char(buffer2) // "'")
    end select

    !> Calculate relaxed density of SSR or SA-REKS state
    call getChildValue(node, "RelaxedDensity", ctrl%reksInp%tRD, default=.false.)
    !> Calculate nonadiabatic coupling vectors
    call getChildValue(node, "NonAdiabaticCoupling", ctrl%reksInp%tNAC, default=.false.)

    !> Print level in standard output file
    call getChildValue(node, "VerbosityLevel", ctrl%reksInp%Plevel, default=1)

  end subroutine readSSR22


  !> Reads SpinTuning block in REKS input
  subroutine readSpinTuning(node, ctrl, nType)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Number of types for atoms
    integer, intent(in) :: nType

    type(fnode), pointer :: value1, child
    type(string) :: buffer, modifier
    type(TListRealR1) :: realBuffer
    integer :: nAtom, iType
    real(dp), allocatable :: tmpTuning(:,:)

    call getChildValue(node, "SpinTuning", value1, "", child=child, &
        & modifier=modifier, allowEmptyValue=.true.)
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
      call getChildValue(child, "", 1, realBuffer, modifier=modifier)
      nAtom = len(realBuffer)
      if (nAtom /= nType) then
        call detailedError(node, "Incorrect number of 'SpinTuning' block: " &
            & // i2c(nAtom) // " supplied, " &
            & // i2c(nType) // " required.")
      end if
      allocate(tmpTuning(1,nAtom))
      call asArray(realBuffer, tmpTuning)
      call destruct(realBuffer)
      allocate(ctrl%reksInp%Tuning(nType))
      ctrl%reksInp%Tuning(:) = tmpTuning(1,:)
    end if

  end subroutine readSpinTuning


  !> Parses Chimes related options.
  subroutine parseChimes(root, chimesRepInput)
    type(fnode), pointer, intent(in) :: root
    type(TChimesRepInp), allocatable, intent(out) :: chimesRepInput

    type(fnode), pointer :: chimes
    type(string) :: buffer

  #:if WITH_CHIMES
    type(string), allocatable :: searchPath(:)
  #:endif

    character(len=:), allocatable :: chimesFile

    call getChild(root, "Chimes", chimes, requested=.false.)
    if (.not. associated(chimes)) return
    #:if WITH_CHIMES
      allocate(chimesRepInput)
      call getChildValue(chimes, "ParameterFile", buffer, default="chimes.dat")
      chimesFile = unquote(char(buffer))
      call getParamSearchPath(searchPath)
      call findFile(searchPath, chimesFile, chimesRepInput%chimesFile)
      if (.not. allocated(chimesRepInput%chimesFile)) then
        call error("Could not find ChIMES parameter file '" // chimesFile // "'")
      end if
    #:else
      call detailedError(chimes, "ChIMES repuslive correction requested, but code was compiled&
          & without ChIMES support")
    #:endif

  end subroutine parseChimes


  !> Returns parser version for a given input version or throws an error if not possible.
  function parserVersionFromInputVersion(versionString, node) result(parserVersion)

    !> Input version string
    character(len=*), intent(in) :: versionString

    !> Input version node (needed for error messagess)
    type(fnode), pointer :: node

    !> Corresponding parser version.
    integer :: parserVersion

    integer :: ii

    do ii = 1, size(versionMaps)
      if (versionMaps(ii)%inputVersion == versionString) then
        parserVersion = versionMaps(ii)%parserVersion
        return
      end if
    end do

    call detailedError(node, "Program version '"// trim(versionString) // "' is not recognized")

  end function parserVersionFromInputVersion


end module dftbp_dftbplus_parser
