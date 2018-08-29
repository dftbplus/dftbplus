!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Fills the derived type with the input parameters from an HSD or an XML file.
module parser
  use globalenv
  use assert
  use accuracy
  use constants
  use inputdata_module
  use typegeometryhsd
  use hsdparser, only : dumpHSD, dumpHSDAsXML, getNodeHSDName
  use hsdutils
  use hsdutils2
  use charmanip
  use message
  use linkedlist
  use unitconversion
  use oldcompat
  use lapackroutines, only : matinv
  use periodic
  use simplealgebra, only: determinant33
  use dispersions
  use dftbplusu
  use slakocont
  use slakoeqgrid
  use repcont
  use repspline
  use reppoly
  use commontypes
  use oldskdata
  use xmlf90
#:if WITH_SOCKETS
  use ipisocket, only : IPI_PROTOCOLS
#:endif
  implicit none

  private

  public :: parseHsdInput, parserVersion


  ! Default file names

  !> Main HSD input file
  character(len=*), parameter :: hsdInputName = "dftb_in.hsd"

  !> XML input file
  character(len=*), parameter :: xmlInputName = "dftb_in.xml"

  !> Processed HSD input
  character(len=*), parameter :: hsdProcInputName = "dftb_pin.hsd"

  !> Processed  XML input
  character(len=*), parameter :: xmlProcInputName = "dftb_pin.xml"

  !> Tag at the head of the input document tree
  character(len=*), parameter :: rootTag = "dftb_in"


  !> Version of the current parser
  integer, parameter :: parserVersion = 6


  !> Version of the oldest parser for which compatibility is still maintained
  integer, parameter :: minVersion = 1


  !> Container type for parser related flags.
  type TParserFlags

    !> stop after parsing?
    logical :: tStop

    !> Continue despite unprocessed nodes
    logical :: tIgnoreUnprocessed

    !> XML output?
    logical :: tWriteXML

    !> HSD output?
    logical :: tWriteHSD
  end type TParserFlags

contains


  !> Parse input from an HSD/XML file
  subroutine parseHsdInput(input)

    !> Returns initialised input variables on exit
    type(inputData), intent(out) :: input

    type(fnode), pointer :: hsdTree
    type(fnode), pointer :: root, tmp, hamNode, child, dummy
    type(TParserflags) :: parserFlags
    logical :: tHSD, missing

    write(stdOut, "(/, A, /)") "***  Parsing and initializing"

    ! Read in the input
    call readHSDOrXML(hsdInputName, xmlInputName, rootTag, hsdTree, tHSD, &
        &missing)

    !! If input is missing return
    if (missing) then
      call error("No input file found.")
    end if

    write(stdout, '(A,1X,I0,/)') 'Parser version:', parserVersion
    if (tHSD) then
      write(stdout, "(A)") "Interpreting input file '" // hsdInputName // "'"
    else
      write(stdout, "(A)") "Interpreting input file '" // xmlInputName //  "'"
    end if
    write(stdout, "(A)") repeat("-", 80)

    ! Get the root of all evil ;-)
    call getChild(hsdTree, rootTag, root)

    ! Handle parser options
    call getChildValue(root, "ParserOptions", dummy, "", child=child, &
        &list=.true., allowEmptyValue=.true., dummyValue=.true.)
    call readParserOptions(child, root, parserFlags)

    ! Read in the different blocks

    ! Atomic geometry and boundary conditions
    call getChild(root, "Geometry", tmp)
    call readGeometry(tmp, input)

    ! electronic Hamiltonian
    call getChildValue(root, "Hamiltonian", hamNode)
    call readHamiltonian(hamNode, input%ctrl, input%geom, input%slako)

    ! Geometry driver
    call getChildValue(root, "Driver", tmp, "", child=child, allowEmptyValue=.true.)
    call readDriver(tmp, child, input%geom, input%ctrl)

    ! excited state options
    call getChildValue(root, "ExcitedState", dummy, "", child=child, list=.true., &
        & allowEmptyValue=.true., dummyValue=.true.)
    call readExcited(child, input%ctrl)

    ! Analysis of properties
    call getChildValue(root, "Analysis", dummy, "", child=child, list=.true., &
        & allowEmptyValue=.true., dummyValue=.true.)
    call readAnalysis(child, input%ctrl, input%geom)

    ! Options for calculation
    call getChildValue(root, "Options", dummy, "", child=child, list=.true., &
        & allowEmptyValue=.true., dummyValue=.true.)
    call readOptions(child, input%ctrl)

    ! Read W values if needed by Hamitonian or excited state calculation
    call readSpinConstants(hamNode, input%geom, input%slako, input%ctrl)

    call readParallel(root, input%ctrl%parallelOpts)

    ! input data strucutre has been initialised
    input%tInitialized = .true.

    ! Issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, parserFlags%tIgnoreUnprocessed)

    ! Dump processed tree in HSD and XML format
    if (tIoProc .and. parserFlags%tWriteHSD) then
      call dumpHSD(hsdTree, hsdProcInputName)
      write(stdout, '(/,/,A)') "Processed input in HSD format written to '" &
          &// hsdProcInputName // "'"
    end if
    if (tIoProc .and. parserFlags%tWriteXML) then
      call dumpHSDAsXML(hsdTree, xmlProcInputName)
      write(stdout, '(A,/)') "Processed input in XML format written to '" &
          &// xmlProcInputName // "'"
    end if

    ! Stop, if only parsing is required
    if (parserFlags%tStop) then
      call error("Keyword 'StopAfterParsing' is set to Yes. Stopping.")
    end if

    call destroyNode(hsdTree)

  end subroutine parseHsdInput


  !> Read in parser options (options not passed to the main code)
  subroutine readParserOptions(node, root, flags)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Root of the entire tree (in case it needs to be converted, for example because of compability
    !> options)
    type(fnode), pointer :: root

    !> Contains parser flags on exit.
    type(TParserFlags), intent(out) :: flags

    integer :: inputVersion
    type(fnode), pointer :: child

    ! Check if input needs compatibility conversion.
    call getChildValue(node, "ParserVersion", inputVersion, parserVersion, &
        &child=child)
    if (inputVersion < 1 .or. inputVersion > parserVersion) then
      call detailedError(child, "Invalid parser version (" // i2c(inputVersion)&
          &// ")")
    elseif (inputVersion < minVersion) then
      call detailedError(child, &
          &"Sorry, no compatibility mode for parser version " &
          &// i2c(inputVersion) // " (too old)")
    elseif (inputVersion /= parserVersion) then
      write(stdout, "(A,I2,A,I2,A)") "***  Converting input from version ", &
          &inputVersion, " to version ", parserVersion, " ..."
      call convertOldHSD(root, inputVersion, parserVersion)
      write(stdout, "(A,/)") "***  Done."
    end if

    call getChildValue(node, "WriteHSDInput", flags%tWriteHSD, .true.)
    call getChildValue(node, "WriteXMLInput", flags%tWriteXML, .false.)
    if (.not. (flags%tWriteHSD .or. flags%tWriteXML)) then
      call detailedWarning(node, &
          &"WriteHSDInput and WriteXMLInput both turned off. You are not&
          & guaranteed" &
          &// newline // &
          &" to able to obtain the same results with a later version of the&
          & code!")
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
    type(inputData), intent(inout) :: input

    type(fnode), pointer :: value, child
    type(string) :: buffer

    call getChildValue(node, "", value, child=child)
    call getNodeName(value, buffer)
    select case (char(buffer))
    case ("genformat")
      call readTGeometryGen(value, input%geom)
    case default
      call setUnprocessed(value)
      call readTGeometryHSD(child, input%geom)
    end select

  end subroutine readGeometry


  !> Read in driver properties
  subroutine readDriver(node, parent, geom, ctrl)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Parent of node (for error messages)
    type(fnode), pointer :: parent

    !> Control structure to be filled
    type(TGeometry), intent(in) :: geom

    !> Nr. of atoms in the system
    type(control), intent(inout) :: ctrl

    type(fnode), pointer :: child, child2, child3, value, value2, field

    type(string) :: buffer, buffer2, modifier
#:if WITH_SOCKETS
    character(lc) :: sTmp
#:endif

    ctrl%tGeoOpt = .false.
    ctrl%tCoordOpt = .false.
    ctrl%tLatOpt = .false.

    ctrl%iGeoOpt = 0
    ctrl%tMD = .false.
    ctrl%iThermostat = 0
    ctrl%tForces = .false.
    ctrl%tSetFillingTemp = .false.

    call getNodeName2(node, buffer)
    driver: select case (char(buffer))
    case ("")
      continue
    case ("none")
      continue
    case ("steepestdescent")
      ! Steepest downhill optimisation

      ctrl%iGeoOpt = 1
      ctrl%tForces = .true.
      ctrl%restartFreq = 1

      call getChildValue(node, "LatticeOpt", ctrl%tLatOpt, .false.)
      if (ctrl%tLatOpt) then
        call getChildValue(node, "Pressure", ctrl%pressure, 0.0_dp, &
            & modifier=modifier, child=child)
        call convertByMul(char(modifier), pressureUnits, child, &
            & ctrl%pressure)
        call getChildValue(node, "FixAngles", ctrl%tLatOptFixAng, .false.)
        if (ctrl%tLatOptFixAng) then
          call getChildValue(node, "FixLengths", ctrl%tLatOptFixLen, &
              & (/.false.,.false.,.false./))
        else
          call getChildValue(node, "Isotropic", ctrl%tLatOptIsotropic, .false.)
        end if
        call getChildValue(node, "MaxLatticeStep", ctrl%maxLatDisp, 0.2_dp)
      end if
      call getChildValue(node, "MovedAtoms", buffer2, "1:-1", child=child, &
          &multiple=.true.)
      call convAtomRangeToInt(char(buffer2), geom%speciesNames, geom%species, &
          &child, ctrl%indMovedAtom)

      ctrl%nrMoved = size(ctrl%indMovedAtom)
      ctrl%tCoordOpt = (ctrl%nrMoved /= 0)
      if (ctrl%tCoordOpt) then
        call getChildValue(node, "MaxAtomStep", ctrl%maxAtomDisp, 0.2_dp)
      end if
      call getChildValue(node, "MaxForceComponent", ctrl%maxForce, 1e-4_dp, &
          &modifier=modifier, child=field)
      call convertByMul(char(modifier), forceUnits, field, ctrl%maxForce)
      call getChildValue(node, "MaxSteps", ctrl%maxRun, 200)
      call getChildValue(node, "StepSize", ctrl%deltaT, 100.0_dp, &
          &modifier=modifier, child=field)
      call convertByMul(char(modifier), timeUnits, field, ctrl%deltaT)
      call getChildValue(node, "OutputPrefix", buffer2, "geo_end")
      ctrl%outFile = unquote(char(buffer2))
      call getChildValue(node, "AppendGeometries", ctrl%tAppendGeo, .false.)
      call getChildValue(node, "ConvergentForcesOnly", ctrl%tConvrgForces, &
          & .true.)
      call readGeoConstraints(node, ctrl, geom%nAtom)
      if (ctrl%tLatOpt) then
        if (ctrl%nrConstr/=0) then
          call error("Lattice optimisation and constraints currently&
              & incompatible.")
        end if
        if (ctrl%nrMoved/=0.and.ctrl%nrMoved<geom%nAtom) then
          call error("Subset of optimising atoms not currently possible with&
              & lattice optimisation.")
        end if
      end if
      ctrl%tGeoOpt = ctrl%tLatOpt .or. ctrl%tCoordOpt

    case ("conjugategradient")
      ! Conjugate gradient location optimisation

      ctrl%iGeoOpt = 2
      ctrl%tForces = .true.
      ctrl%restartFreq = 1
      call getChildValue(node, "LatticeOpt", ctrl%tLatOpt, .false.)
      if (ctrl%tLatOpt) then
        call getChildValue(node, "Pressure", ctrl%pressure, 0.0_dp, &
            & modifier=modifier, child=child)
        call convertByMul(char(modifier), pressureUnits, child, &
            & ctrl%pressure)
        call getChildValue(node, "FixAngles", ctrl%tLatOptFixAng, .false.)
        if (ctrl%tLatOptFixAng) then
          call getChildValue(node, "FixLengths", ctrl%tLatOptFixLen, &
              & (/.false.,.false.,.false./))
        else
          call getChildValue(node, "Isotropic", ctrl%tLatOptIsotropic, .false.)
        end if
        call getChildValue(node, "MaxLatticeStep", ctrl%maxLatDisp, 0.2_dp)
      end if
      call getChildValue(node, "MovedAtoms", buffer2, "1:-1", child=child, &
          &multiple=.true.)
      call convAtomRangeToInt(char(buffer2), geom%speciesNames, geom%species, &
          &child, ctrl%indMovedAtom)

      ctrl%nrMoved = size(ctrl%indMovedAtom)
      ctrl%tCoordOpt = (ctrl%nrMoved /= 0)
      if (ctrl%tCoordOpt) then
        call getChildValue(node, "MaxAtomStep", ctrl%maxAtomDisp, 0.2_dp)
      end if
      call getChildValue(node, "MaxForceComponent", ctrl%maxForce, 1e-4_dp, &
          &modifier=modifier, child=field)
      call convertByMul(char(modifier), forceUnits, field, ctrl%maxForce)
      call getChildValue(node, "MaxSteps", ctrl%maxRun, 200)
      call getChildValue(node, "OutputPrefix", buffer2, "geo_end")
      ctrl%outFile = unquote(char(buffer2))
      call getChildValue(node, "AppendGeometries", ctrl%tAppendGeo, .false.)
      call getChildValue(node, "ConvergentForcesOnly", ctrl%tConvrgForces, &
          & .true.)
      call readGeoConstraints(node, ctrl, geom%nAtom)
      if (ctrl%tLatOpt) then
        if (ctrl%nrConstr/=0) then
          call error("Lattice optimisation and constraints currently&
              & incompatible.")
        end if
        if (ctrl%nrMoved/=0.and.ctrl%nrMoved<geom%nAtom) then
          call error("Subset of optimising atoms not currently possible with&
              & lattice optimisation.")
        end if
      end if
      ctrl%tGeoOpt = ctrl%tLatOpt .or. ctrl%tCoordOpt

    case("gdiis")
      ! Gradient DIIS optimisation, only stable in the quadratic region

      ctrl%iGeoOpt = 3
      ctrl%tForces = .true.
      ctrl%restartFreq = 1
      call getChildValue(node, "alpha", ctrl%deltaGeoOpt, 1.0E-1_dp)
      call getChildValue(node, "Generations", ctrl%iGenGeoOpt, 8)
      call getChildValue(node, "LatticeOpt", ctrl%tLatOpt, .false.)
      if (ctrl%tLatOpt) then
        call getChildValue(node, "Pressure", ctrl%pressure, 0.0_dp, &
            & modifier=modifier, child=child)
        call convertByMul(char(modifier), pressureUnits, child, &
            & ctrl%pressure)
        call getChildValue(node, "FixAngles", ctrl%tLatOptFixAng, .false.)
        if (ctrl%tLatOptFixAng) then
          call getChildValue(node, "FixLengths", ctrl%tLatOptFixLen, &
              & (/.false.,.false.,.false./))
        else
          call getChildValue(node, "Isotropic", ctrl%tLatOptIsotropic, .false.)
        end if
        call getChildValue(node, "MaxLatticeStep", ctrl%maxLatDisp, 0.2_dp)
      end if
      call getChildValue(node, "MovedAtoms", buffer2, "1:-1", child=child, &
          &multiple=.true.)
      call convAtomRangeToInt(char(buffer2), geom%speciesNames, geom%species, &
          &child, ctrl%indMovedAtom)

      ctrl%nrMoved = size(ctrl%indMovedAtom)
      ctrl%tCoordOpt = (ctrl%nrMoved /= 0)
      call getChildValue(node, "MaxForceComponent", ctrl%maxForce, 1e-4_dp, &
          &modifier=modifier, child=field)
      call convertByMul(char(modifier), forceUnits, field, ctrl%maxForce)
      call getChildValue(node, "MaxSteps", ctrl%maxRun, 200)
      call getChildValue(node, "OutputPrefix", buffer2, "geo_end")
      ctrl%outFile = unquote(char(buffer2))
      call getChildValue(node, "AppendGeometries", ctrl%tAppendGeo, .false.)
      call getChildValue(node, "ConvergentForcesOnly", ctrl%tConvrgForces, &
          & .true.)
      call readGeoConstraints(node, ctrl, geom%nAtom)
      if (ctrl%tLatOpt) then
        if (ctrl%nrConstr/=0) then
          call error("Lattice optimisation and constraints currently&
              & incompatible.")
        end if
        if (ctrl%nrMoved/=0.and.ctrl%nrMoved<geom%nAtom) then
          call error("Subset of optimising atoms not currently possible with&
              & lattice optimisation.")
        end if
      end if
      ctrl%tGeoOpt = ctrl%tLatOpt .or. ctrl%tCoordOpt

    case ("lbfgs")

      ctrl%iGeoOpt = 4

      ctrl%tForces = .true.
      ctrl%restartFreq = 1
      call getChildValue(node, "LatticeOpt", ctrl%tLatOpt, .false.)
      if (ctrl%tLatOpt) then
        call getChildValue(node, "Pressure", ctrl%pressure, 0.0_dp, &
            & modifier=modifier, child=child)
        call convertByMul(char(modifier), pressureUnits, child, &
            & ctrl%pressure)
        call getChildValue(node, "FixAngles", ctrl%tLatOptFixAng, .false.)
        if (ctrl%tLatOptFixAng) then
          call getChildValue(node, "FixLengths", ctrl%tLatOptFixLen, &
              & (/.false.,.false.,.false./))
        else
          call getChildValue(node, "Isotropic", ctrl%tLatOptIsotropic, .false.)
        end if
        call getChildValue(node, "MaxLatticeStep", ctrl%maxLatDisp, 0.2_dp)
      end if
      call getChildValue(node, "MovedAtoms", buffer2, "1:-1", child=child, &
          &multiple=.true.)
      call convAtomRangeToInt(char(buffer2), geom%speciesNames, geom%species, &
          &child, ctrl%indMovedAtom)

      ctrl%nrMoved = size(ctrl%indMovedAtom)
      ctrl%tCoordOpt = (ctrl%nrMoved /= 0)
      if (ctrl%tCoordOpt) then
        call getChildValue(node, "MaxAtomStep", ctrl%maxAtomDisp, 0.2_dp)
      end if
      call getChildValue(node, "MaxForceComponent", ctrl%maxForce, 1e-4_dp, &
          &modifier=modifier, child=field)
      call convertByMul(char(modifier), forceUnits, field, ctrl%maxForce)
      call getChildValue(node, "MaxSteps", ctrl%maxRun, 200)
      call getChildValue(node, "OutputPrefix", buffer2, "geo_end")
      ctrl%outFile = unquote(char(buffer2))
      call getChildValue(node, "AppendGeometries", ctrl%tAppendGeo, .false.)
      call getChildValue(node, "ConvergentForcesOnly", ctrl%tConvrgForces, &
          & .true.)
      call readGeoConstraints(node, ctrl, geom%nAtom)
      if (ctrl%tLatOpt) then
        if (ctrl%nrConstr/=0) then
          call error("Lattice optimisation and constraints currently&
              & incompatible.")
        end if
        if (ctrl%nrMoved/=0.and.ctrl%nrMoved<geom%nAtom) then
          call error("Subset of optimising atoms not currently possible with&
              & lattice optimisation.")
        end if
      end if
      ctrl%tGeoOpt = ctrl%tLatOpt .or. ctrl%tCoordOpt

      allocate(ctrl%lbfgsInp)
      call getChildValue(node, "Memory", ctrl%lbfgsInp%memory, 20)

    case("secondderivatives")
      ! currently only numerical derivatives of forces is implemented

      ctrl%tDerivs = .true.
      ctrl%tForces = .true.
      call getChildValue(node, "Atoms", buffer2, "1:-1", child=child, &
          &multiple=.true.)
      call convAtomRangeToInt(char(buffer2), geom%speciesNames, geom%species, &
          &child, ctrl%indMovedAtom)
      ctrl%nrMoved = size(ctrl%indMovedAtom)
      if (ctrl%nrMoved == 0) then
        call error("No atoms specified for derivatives calculation.")
      end if
      call getChildValue(node, "Delta", ctrl%deriv2ndDelta, 1.0E-4_dp, &
          & modifier=modifier, child=field)
      call convertByMul(char(modifier), lengthUnits, field, ctrl%deriv2ndDelta)
      ctrl%tConvrgForces = .true.

    case ("velocityverlet")
      ! molecular dynamics

      ctrl%tForces = .true.
      ctrl%tMD = .true.

      call getChildValue(node, "MDRestartFrequency", ctrl%restartFreq, 1)
      call getChildValue(node, "MovedAtoms", buffer2, "1:-1", child=child, &
          &multiple=.true.)
      call convAtomRangeToInt(char(buffer2), geom%speciesNames, geom%species,&
          &child, ctrl%indMovedAtom)
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
      call convertByMul(char(modifier), timeUnits, field, ctrl%deltaT)

      call getChildValue(node, "Thermostat", value, child=child)
      call getNodeName(value, buffer2)

      call getChildValue(node, "ConvergentForcesOnly", ctrl%tConvrgForces, &
          & .true.)

      thermostat: select case(char(buffer2))
      case ("berendsen")
        ctrl%iThermostat = 2
        ! Read temperature or temperature profiles
        call getChildValue(value, "Temperature", value2, modifier=modifier, &
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

        !call getChildValue(value, "CouplingStrength", ctrl%wvScale)
        call getChild(value, "CouplingStrength", child=child2, &
            & requested=.false.)
        if (associated(child2)) then
          call getChildValue(child2, "", ctrl%wvScale)
          call getChild(value, "Timescale",child=child2,modifier=modifier,&
              &requested=.false.)
          if (associated(child2)) call error("Only Coupling strength OR &
              &Timescale can be set for Berendsen thermostats.")
        else
          call getChild(value, "Timescale",child=child2,modifier=modifier,&
              &requested=.false.)
          if (associated(child2)) then
            call getChildValue(child2, "", ctrl%wvScale, &
                & modifier=modifier, child=child3)
            call convertByMul(char(modifier), timeUnits, child3, &
                & ctrl%wvScale)
            ctrl%wvScale = ctrl%deltaT / ctrl%wvScale
          else
            call error("Either CouplingStrength or Timescale must be set&
                & for Berendsen thermostats.")
          end if
        end if

        call getChildValue(value, "AdaptFillingTemp", ctrl%tSetFillingTemp, &
            &.false.)

      case ("nosehoover")
        ctrl%iThermostat = 3
        ! Read temperature or temperature profiles
        call getChildValue(value, "Temperature", value2, modifier=modifier, &
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

        call getChildValue(value, "CouplingStrength", ctrl%wvScale, &
            & modifier=modifier, child=field)
        call convertByMul(char(modifier), freqUnits, field, ctrl%wvScale)

        call getChildValue(value, "ChainLength", ctrl%nh_npart, 3)
        call getChildValue(value, "Order", ctrl%nh_nys, 3)
        call getChildValue(value, "IntegratorSteps", ctrl%nh_nc, 1)

        call getChild(value, "Restart",  child=child3, requested=.false.)
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

        call getChildValue(value, "AdaptFillingTemp", ctrl%tSetFillingTemp, &
            &.false.)

      case ("andersen")
        ctrl%iThermostat = 1
        ! Read temperature or temperature profiles
        call getChildValue(value, "Temperature", value2, modifier=modifier, &
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

        call getChildValue(value, "ReselectProbability", ctrl%wvScale, &
            &child=child3)
        if (ctrl%wvScale <= 0.0_dp .or. ctrl%wvScale > 1.0_dp) then
          call detailedError(child3, &
              &"ReselectProbability must be in the range (0,1]!")
        end if
        call getChildValue(value, "ReselectIndividually", ctrl%tRescale)
        call getChildValue(value, "AdaptFillingTemp", ctrl%tSetFillingTemp, &
            &.false.)

      case ("none")
        ctrl%iThermostat = 0
        allocate(ctrl%tempSteps(1))
        allocate(ctrl%tempValues(1))

        if (ctrl%tReadMDVelocities) then
          ! without a thermostat, if we know the initial velocities, we do not
          ! need a temperature, so just set it to something 'safe'
          ctrl%tempAtom = minTemp
        else
          call getChildValue(value, "InitialTemperature", ctrl%tempAtom, &
              &modifier=modifier, child=field)
          if (ctrl%tempAtom < 0.0_dp) then
            call detailedError(field, "Negative temperature")
          end if
          call convertByMul(char(modifier), energyUnits, field, ctrl%tempAtom)
          if (ctrl%tempAtom < minTemp) then
            ctrl%tempAtom = minTemp
          end if
        end if
      case default
        call getNodeHSDName(value, buffer2)
        call detailedError(child, "Invalid thermostat '" // char(buffer2) // "'")
      end select thermostat

      if (ctrl%maxRun < -1) then
        call getChildValue(node, "Steps", ctrl%maxRun)
      end if

      call getChildValue(node, "OutputPrefix", buffer2, "geo_end")
      ctrl%outFile = unquote(char(buffer2))

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
          call convertByMul(char(modifier), pressureUnits, child2, &
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
              call convertByMul(char(modifier), timeUnits, child3, &
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

      call readXlbomdOptions(node, ctrl%xlbomd)

      call getInputMasses(node, geom, ctrl%masses)

    case ("socket")
      ! external socket control of the run (once initialised from input)
#:if WITH_SOCKETS
      ctrl%tForces = .true.
      allocate(ctrl%socketInput)
      call getChild(node, 'File', child=child2, requested=.false.)
      call getChild(node, 'Host', child=child3, requested=.false.)
      if (associated(child2) .eqv. associated(child3)) then
        call error('Either Host or File (but not both) must be set for socket&
            & communication')
      end if

      ! File communiation
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

      call getChildValue(node, "Protocol", value, "i-PI", child=child)
      call getNodeName(value, buffer)
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

  end subroutine readDriver


  !> Extended lagrangian options
  subroutine readXlbomdOptions(node, input)

    !> node in the input tree
    type(fnode), pointer :: node

    !> extracted settings on exit
    type(XlbomdInp), allocatable, intent(out) :: input

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
      call getChildValue(pRoot, 'SccTolerance', input%sccTol, 1e-5_dp)
      input%scale = 1.0_dp
      input%useInverseJacobian = .false.
      input%readInverseJacobian = .false.
    end if

  end subroutine readXlbomdOptions


  !> Reads geometry constraints.
  subroutine readGeoConstraints(node, ctrl, nAtom)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(control), intent(inout) :: ctrl

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    type(fnode), pointer :: value, child
    type(string) :: buffer
    type(listIntR1) :: intBuffer
    type(listRealR1) :: realBuffer

    call getChildValue(node, "Constraints", value, "", child=child, &
        &allowEmptyValue=.true.)
    call getNodeName2(value, buffer)
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
    type(control), intent(inout) :: ctrl

    !> Total number of all atoms
    integer, intent(in) :: nAtom

    type(fnode), pointer :: value, child
    type(string) :: buffer, modifier
    type(listRealR1) :: realBuffer
    integer :: nVelocities
    real(dp), allocatable :: tmpVelocities(:,:)

    call getChildValue(node, "Velocities", value, "", child=child, &
        & modifier=modifier, allowEmptyValue=.true.)
    call getNodeName2(value, buffer)
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
        call convertByMul(char(modifier), VelocityUnits, child, &
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
      return
    end if

    allocate(masses(geo%nAtom))
    masses(:) = -1.0_dp
    do ii = 1, getLength(children)
      call getItem1(children, ii, child2)
      call getChildValue(child2, "Atoms", buffer, child=child3, multiple=.true.)
      call convAtomRangeToInt(char(buffer), geo%speciesNames, geo%species, child3, pTmpI1)
      call getChildValue(child2, "MassPerAtom", rTmp, modifier=modifier, child=child)
      call convertByMul(char(modifier), massUnits, child, rTmp)
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
  subroutine readHamiltonian(node, ctrl, geo, slako)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(control), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Slater-Koster structure to be filled
    type(slater), intent(inout) :: slako

    type(string) :: buffer

    call getNodeName(node, buffer)
    select case (char(buffer))
    case ("dftb")
      call readDFTBHam(node, ctrl, geo, slako)
    case default
      call detailedError(node, "Invalid Hamiltonian")
    end select

  end subroutine readHamiltonian


  !> Reads DFTB-Hamiltonian
  subroutine readDFTBHam(node, ctrl, geo, slako)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(control), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Slater-Koster structure to be filled
    type(slater), intent(inout) :: slako

    type(fnode), pointer :: value, value2, child, child2, child3, field
    type(fnodeList), pointer :: children
    type(string) :: buffer, buffer2, modifier
    type(listInt) :: li
    type(listIntR1) :: li1
    type(listRealR1) :: lr1
    type(listInt), allocatable :: liN(:)
    type(listIntR1), allocatable :: li1N(:)
    type(listReal), allocatable :: lrN(:)
    type(listCharLc), allocatable :: skFiles(:,:)
    type(listString) :: lStr
    type(listIntR1), allocatable :: angShells(:)
    type(listRealR2) :: lCharges
    type(listRealR1) :: lBlurs
    logical, allocatable :: repPoly(:,:)
    integer :: iSp1, iSp2, iSh1, ii, jj, kk, ind
    character(lc) :: prefix, suffix, separator, elem1, elem2, strTmp
    character(lc) :: errorStr
    logical :: tLower, tExist
    real(dp), allocatable :: kpts(:,:)
    integer, allocatable :: tmpI1(:)
    integer, allocatable :: pTmpI1(:)
    real(dp), allocatable :: tmpR1(:)
    real(dp), allocatable :: tmpR2(:,:)
    real(dp) :: rTmp, rTmp3(3)
    integer, allocatable :: iTmpN(:)
    real(dp) :: coeffsAndShifts(3, 4)
    integer :: nShell, skInterMeth
    character(1) :: tmpCh
    logical :: tShellIncl(4), tFound
    integer :: angShell(maxL+1), angShellOrdered(maxL+1)
    integer :: fp, iErr
    logical :: tBadIntegratingKPoints
    integer :: nElem

    ! Read in maximal angular momenta or selected shells
    do ii = 1, maxL+1
      angShellOrdered(ii) = ii - 1
    end do
    call getChild(node, "MaxAngularMomentum", child)
    allocate(angShells(geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      call init(angShells(iSp1))
      call getChildValue(child, geo%speciesNames(iSp1), value, child=child2)
      call getNodeName(value, buffer)
      select case(char(buffer))
      case("selectedshells")
        call init(lStr)
        call getChildValue(value, "", lStr)
        do ii = 1, len(lStr)
          call get(lStr, strTmp, ii)
          strTmp = tolower(unquote(trim(strTmp)))
          if (len_trim(strTmp) > 4 .or. len_trim(strTmp) < 1) then
            call detailedError(value, "Invalid shell selection '" &
                &// trim(strTmp) &
                &// "'. Nr. of selected shells must be between 1 and 4.")
          end if
          tShellIncl(:) = .false.
          nShell = len_trim(strTmp)
          do jj = 1, nShell
            tmpCh = strTmp(jj:jj)
            tFound = .false.
            do kk = 1, size(orbitalNames)
              if (tmpCh == trim(orbitalNames(kk))) then
                if (tShellIncl(kk)) then
                  call detailedError(value, "Double selection of the same shell&
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
              call detailedError(value, "Invalid shell name '" // tmpCh // "'")
            end if
          end do
          call append(angShells(iSp1), angShell(1:nShell))
        end do
        call destruct(lStr)

      case(textNodeName)
        call getChildValue(child2, "", buffer)
        strTmp = unquote(char(buffer))
        do jj = 1, size(orbitalNames)
          if (trim(strTmp) == trim(orbitalNames(jj))) then
            call append(angShells(iSp1), angShellOrdered(:jj))
          end if
        end do
        if (len(angShells(iSp1)) < 1) then
          call detailedError(child2, "Invalid orbital name '" // &
              &trim(strTmp) // "'")
        end if

      case default
        call getNodeHSDName(value, buffer)
        call detailedError(child2, "Invalid shell specification method '" //&
            & char(buffer) // "'")
      end select
    end do

    ! Orbitals and angular momenta for the given shells (once the SK files will contain the full
    ! information about the basis, this will be moved to the SK reading routine).
    allocate(slako%orb)
    allocate(slako%orb%nShell(geo%nSpecies))
    allocate(slako%orb%nOrbSpecies(geo%nSpecies))
    allocate(slako%orb%nOrbAtom(geo%nAtom))
    slako%orb%mOrb = 0
    slako%orb%mShell = 0
    do iSp1 = 1, geo%nSpecies
      slako%orb%nShell(iSp1) = 0
      slako%orb%nOrbSpecies(iSp1) = 0
      do ii = 1, len(angShells(iSp1))
        call intoArray(angShells(iSp1), angShell, nShell, ii)
        slako%orb%nShell(iSp1) = slako%orb%nShell(iSp1) + nShell
        do jj = 1, nShell
          slako%orb%nOrbSpecies(iSp1) = slako%orb%nOrbSpecies(iSp1) &
              &+ 2 * angShell(jj) + 1
        end do
      end do
    end do
    slako%orb%mShell = maxval(slako%orb%nShell)
    slako%orb%mOrb = maxval(slako%orb%nOrbSpecies)
    slako%orb%nOrbAtom(:) = slako%orb%nOrbSpecies(geo%species(:))
    slako%orb%nOrb = sum(slako%orb%nOrbAtom)

    allocate(slako%orb%angShell(slako%orb%mShell, geo%nSpecies))
    allocate(slako%orb%iShellOrb(slako%orb%mOrb, geo%nSpecies))
    allocate(slako%orb%posShell(slako%orb%mShell+1, geo%nSpecies))
    slako%orb%angShell(:,:) = 0
    do iSp1 = 1, geo%nSpecies
      ind = 1
      iSh1 = 1
      do ii = 1, len(angShells(iSp1))
        call intoArray(angShells(iSp1), angShell, nShell, ii)
        do jj = 1, nShell
          slako%orb%posShell(iSh1, iSp1) = ind
          slako%orb%angShell(iSh1, iSp1) = angShell(jj)
          slako%orb%iShellOrb(ind:ind+2*angShell(jj), iSp1) = iSh1
          ind = ind + 2 * angShell(jj) + 1
          iSh1 = iSh1 + 1
        end do
        slako%orb%posShell(iSh1, iSp1) = ind
      end do
    end do

    ! Slater-Koster files
    allocate(skFiles(geo%nSpecies, geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
      do iSp2 = 1, geo%nSpecies
        call init(skFiles(iSp2, iSp1))
      end do
    end do
    call getChildValue(node, "SlaterKosterFiles", value, child=child)
    call getNodeName(value, buffer)
    select case(char(buffer))
    case ("type2filenames")
      call getChildValue(value, "Prefix", buffer2, "")
      prefix = unquote(char(buffer2))
      call getChildValue(value, "Suffix", buffer2, "")
      suffix = unquote(char(buffer2))
      call getChildValue(value, "Separator", buffer2, "")
      separator = unquote(char(buffer2))
      call getChildValue(value, "LowerCaseTypeName", tLower, .false.)
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
          call append(skFiles(iSp2, iSp1), strTmp)
          inquire(file=strTmp, exist=tExist)
          if (.not. tExist) then
            call detailedError(value, "SK file with generated name '" &
                &// trim(strTmp) // "' does not exist.")
          end if
        end do
      end do
    case default
      call setUnprocessed(value)
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
    call getChildValue(node, "PolynomialRepulsive", value, "", child=child, &
        &list=.true., allowEmptyValue=.true., dummyValue=.true.)
    call getNodeName2(value, buffer)
    select case (char(buffer))
    case ("")
      repPoly(:,:) = .false.
    case("setforall")
      call getChildValue(value, "", repPoly(1,1))
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
        call detailedError(value, "Assymetric definition (both A-B and B-A must&
            & be defined for using polynomial repulsive)")
      end if
    end select

    call getChildValue(node, "OrbitalResolvedSCC", ctrl%tOrbResolved, .false.)
    call getChildValue(node, "OldSKInterpolation", ctrl%oldSKInter, .false.)
    if (ctrl%oldSKInter) then
      skInterMeth = skEqGridOld
    else
      skInterMeth = skEqGridNew
    end if
    call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, &
        &angShells, ctrl%tOrbResolved, skInterMeth, repPoly)

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
    call getChildValue(node, "SCC", ctrl%tSCC, .false.)
    ifSCC: if (ctrl%tSCC) then
      call getChildValue(node, "ReadInitialCharges", ctrl%tReadChrg, .false.)
      if (.not. ctrl%tReadChrg) then
        call getInitialCharges(node, geo, ctrl%initialCharges)
      end if
      call getChildValue(node, "SCCTolerance", ctrl%sccTol, 1.0e-5_dp)
      call getChildValue(node, "Mixer", value, "Broyden", child=child)
      call getNodeName(value, buffer)
      select case(char(buffer))

      case ("broyden")
        ctrl%iMixSwitch = 3
        call getChildValue(value, "MixingParameter", ctrl%almix, 0.2_dp)
        call getChildValue(value, "InverseJacobiWeight", ctrl%broydenOmega0, &
            &0.01_dp)
        call getChildValue(value, "MinimalWeight", ctrl%broydenMinWeight, &
            &1.0_dp)
        call getChildValue(value, "MaximalWeight", ctrl%broydenMaxWeight, &
            &1.0e5_dp)
        call getChildValue(value, "WeightFactor", ctrl%broydenWeightFac, &
            &1.0e-2_dp)

      case ("anderson")
        ctrl%iMixSwitch = 2
        call getChildValue(value, "MixingParameter", ctrl%almix, 0.05_dp)
        call getChildValue(value, "Generations", ctrl%iGenerations, 4)
        call getChildValue(value, "InitMixingParameter", &
            &ctrl%andersonInitMixing, 0.01_dp)
        call getChildValue(value, "DynMixingParameters", value2, "", &
            &child=child, allowEmptyValue=.true.)
        call getNodeName2(value2, buffer2)
        if (char(buffer2) == "") then
          ctrl%andersonNrDynMix = 0
        else
          call init(lr1)
          call getChildValue(child, "", 2, lr1, child=child2)
          if (len(lr1) < 1) then
            call detailedError(child2, "At least one dynamic mixing parameter&
                & must be defined.")
          end if
          ctrl%andersonNrDynMix = len(lr1)
          allocate(ctrl%andersonDynMixParams(2, ctrl%andersonNrDynMix))
          call asArray(lr1, ctrl%andersonDynMixParams)
          call destruct(lr1)
        end if
        call getChildValue(value, "DiagonalRescaling", ctrl%andersonOmega0, &
            &1.0e-2_dp)

      case ("simple")
        ctrl%iMixSwitch = 1
        call getChildValue(value, "MixingParameter", ctrl%almix, 0.05_dp)

      case("diis")
        ctrl%iMixSwitch = 4
        call getChildValue(value, "InitMixingParameter", ctrl%almix, 0.2_dp)
        call getChildValue(value, "Generations", ctrl%iGenerations, 6)
        call getChildValue(value, "UseFromStart", ctrl%tFromStart, .true.)
      case default
        call getNodeHSDName(value, buffer)
        call detailedError(child, "Invalid mixer '" // char(buffer) // "'")
      end select

      if (geo%tPeriodic) then
        call getChildValue(node, "EwaldParameter", ctrl%ewaldAlpha, 0.0_dp)
        call getChildValue(node, "EwaldTolerance", ctrl%tolEwald, 1.0e-9_dp)
      end if

      call readHCorrection(node, geo, ctrl)

      ! spin
      call getChildValue(node, "SpinPolarisation", value, "", child=child, &
          &allowEmptyValue=.true.)
      call getNodeName2(value, buffer)
      select case(char(buffer))
      case ("")
        ctrl%tSpin = .false.
        ctrl%t2Component = .false.
        ctrl%nrSpinPol = 0.0_dp

      case ("colinear")
        ctrl%tSpin = .true.
        ctrl%t2Component = .false.
        call getChildValue(value, 'UnpairedElectrons', ctrl%nrSpinPol, 0.0_dp)
        call getChildValue(value, 'RelaxTotalSpin', ctrl%tSpinSharedEf, .false.)
        if (.not. ctrl%tReadChrg) then
          call getInitialSpins(value, geo, 1, ctrl%initialSpins)
        end if

      case ("noncolinear")
        ctrl%tSpin = .true.
        ctrl%t2Component = .true.
        if (.not. ctrl%tReadChrg) then
          call getInitialSpins(value, geo, 3, ctrl%initialSpins)
        end if

      case default
        call getNodeHSDName(value, buffer)
        call detailedError(child, "Invalid spin polarisation type '" //&
            & char(buffer) // "'")
      end select

      ctrl%tMulliken = .true.

    end if ifSCC

    ! External electric field
    call getChildValue(node, "ElectricField", value, "", child=child, &
        &allowEmptyValue=.true., dummyValue=.true., list=.true.)

    ! external applied field
    call getChild(child, "External",child2,requested=.false.)
    if (associated(child2)) then
      call getChildValue(child2, "Strength", ctrl%EFieldStrength, &
          & modifier=modifier, child=child3)
      call convertByMul(char(modifier), EFieldUnits, child3, &
          & ctrl%EFieldStrength)
      call getChildValue(child2, "Direction", ctrl%EfieldVector)
      if (sum(ctrl%EfieldVector**2) < 1e-8_dp) then
        call detailedError(child2,"Vector too small")
      else
        ctrl%EfieldVector = ctrl%EfieldVector/sqrt(sum(ctrl%EfieldVector**2))
      end if
      ctrl%tEField = .true.
      ctrl%tMulliken = .true.
      call getChildValue(child2, "Frequency", ctrl%EFieldOmega, 0.0_dp, &
          & modifier=modifier, child=child3)
      call convertByMul(char(modifier), freqUnits, child3, &
          & ctrl%EFieldOmega)
      if (ctrl%EFieldOmega > 0.0) then
        ctrl%EFieldOmega = 2.0_dp * pi * ctrl%EFieldOmega !angular frequency
        ctrl%tTDEfield = .true.
      else
        ctrl%tTDEfield = .false.
        ctrl%EFieldOmega = 0.0_dp
      end if
      ctrl%EfieldPhase = 0
      if (ctrl%tTDEfield) then
        call getChildValue(child2, "Phase", ctrl%EfieldPhase, 0)
      end if
    else
      ctrl%tEField = .false.
    end if

    ! Point charges present
    call getChildren(child, "PointCharges", children)
    if (getLength(children) > 0) then

      if (.not.ctrl%tSCC) then
        call error("External charges can only be used in an SCC calculation")
      end if
      call init(lCharges)
      call init(lBlurs)
      ctrl%nExtChrg = 0
      do ii = 1, getLength(children)
        call getItem1(children, ii, child2)
        call getChildValue(child2, "CoordsAndCharges", value, &
            &modifier=modifier, child=child3)
        call getNodeName(value, buffer)
        select case(char(buffer))
        case (textNodeName)
          call init(lr1)
          call getChildValue(child3, "", 4, lr1, modifier=modifier)
          allocate(tmpR2(4, len(lr1)))
          call asArray(lr1, tmpR2)
          ctrl%nExtChrg = ctrl%nExtChrg + len(lr1)
          call destruct(lr1)
        case ("directread")
          call getChildValue(value, "Records", ind)
          call getChildValue(value, "File", buffer2)
          allocate(tmpR2(4, ind))
          open(newunit=fp, file=unquote(char(buffer2)), form="formatted", status="old",&
              & action="read", iostat=iErr)
          if (iErr /= 0) then
            call detailedError(value, "Could not open file '" &
                &// trim(unquote(char(buffer2))) // "' for direct reading" )
          end if
          read(fp, *, iostat=iErr) tmpR2
          if (iErr /= 0) then
            call detailedError(value, "Error during direct reading '" &
                &// trim(unquote(char(buffer2))) // "'")
          end if
          close(fp)
          ctrl%nExtChrg = ctrl%nExtChrg + ind
        case default
          call detailedError(value, "Invalid block name")
        end select
        call convertByMul(char(modifier), lengthUnits, child3, tmpR2(1:3,:))
        call append(lCharges, tmpR2)
        call getChildValue(child2, "GaussianBlurWidth", rTmp, 0.0_dp, &
            &modifier=modifier, child=child3)
        if (rTmp < 0.0_dp) then
          call detailedError(child3, "Gaussian blur width may not be &
              &negative")
        end if
        call convertByMul(char(modifier), lengthUnits, child3, rTmp)
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
    else
      ctrl%nExtChrg = 0
    end if
    call destroyNodeList(children)

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

      allocate(ctrl%xi(slako%orb%mShell,geo%nSpecies))
      ctrl%xi = 0.0_dp
      do iSp1 = 1, geo%nSpecies
        call getChildValue(child, geo%speciesNames(iSp1), &
            & ctrl%xi(:slako%orb%nShell(iSp1),iSp1), modifier=modifier, &
            & child=child2 )
        call convertByMul(char(modifier), energyUnits, child2, &
            & ctrl%xi(:slako%orb%nShell(iSp1),iSp1))
      end do
    end if

    ! Solver
    call getChildValue(node, "Eigensolver", value, "RelativelyRobust")
    call getNodeName(value, buffer)
    select case(char(buffer))
    case ("qr")
      ctrl%iSolver = 1
    case ("divideandconquer")
      ctrl%iSolver = 2
    case ("relativelyrobust")
      ctrl%iSolver = 3
    end select

    ! Filling (temperature only read, if AdaptFillingTemp was not set for the selected MD
    ! thermostat.)
    call getChildValue(node, "Filling", value, "Fermi", child=child)
    call getNodeName(value, buffer)

    select case (char(buffer))
    case ("fermi")
      ctrl%iDistribFn = 0 ! Fermi function
    case ("methfesselpaxton")
      ! Set the order of the Methfessel-Paxton step function approximation, defaulting to 2
      call getChildValue(value, "Order", ctrl%iDistribFn, 2)
      if (ctrl%iDistribFn < 1) then
        call getNodeHSDName(value, buffer)
        write(errorStr, "(A,A,A,I4)")"Unsuported filling mode for '", &
            & char(buffer),"' :",ctrl%iDistribFn
        call detailedError(child, errorStr)
      end if
    case default
      call getNodeHSDName(value, buffer)
      call detailedError(child, "Invalid filling method '" //char(buffer)// "'")
    end select

    if (.not. ctrl%tSetFillingTemp) then
      call getChildValue(value, "Temperature", ctrl%tempElec, 0.0_dp, &
          &modifier=modifier, child=field)
      call convertByMul(char(modifier), energyUnits, field, ctrl%tempElec)
      if (ctrl%tempElec < minTemp) then
        ctrl%tempElec = minTemp
      end if
    end if

    call getChild(value, "FixedFermiLevel", child=child2, modifier=modifier, &
        & requested=.false.)
    if (associated(child2)) then
      if (ctrl%tSpin .and. .not.ctrl%t2Component) then
        call getChildValue(child2, "", ctrl%Ef(:2), &
            & modifier=modifier, child=child3)
      else
        call getChildValue(child2, "", ctrl%Ef(:1), &
            & modifier=modifier, child=child3)
      end if
      call convertByMul(char(modifier), energyUnits, child3, &
          & ctrl%Ef)
      ctrl%tFixEf = .true.
    else
      ctrl%tFixEf = .false.
    end if

    if (geo%tPeriodic .and. .not.ctrl%tFixEf) then
      call getChildValue(value, "IndependentKFilling", ctrl%tFillKSep, .false.)
    end if

    ! Charge
    call getChildValue(node, "Charge", ctrl%nrChrg, 0.0_dp)

    ! Assume SCC can has usual default number of steps if needed
    tBadIntegratingKPoints = .false.

    ! K-Points
    if (geo%tPeriodic) then
      call getChildValue(node, "KPointsAndWeights", value, child=child, &
          &modifier=modifier)
      call getNodeName(value, buffer)
      select case(char(buffer))

      case ("supercellfolding")
        tBadIntegratingKPoints = .false.
        if (len(modifier) > 0) then
          call detailedError(child, "No modifier is allowed, if the &
              &SupercellFolding scheme is used.")
        end if
        call getChildValue(value, "", coeffsAndShifts)
        if (abs(determinant33(coeffsAndShifts(:,1:3))) - 1.0_dp < -1e-6_dp) then
          call detailedError(value, "Determinant of the supercell matrix must &
              &be greater than 1")
        end if
        if (any(abs(modulo(coeffsAndShifts(:,1:3) + 0.5_dp, 1.0_dp) - 0.5_dp) &
            &> 1e-6_dp)) then
          call detailedError(value, "The components of the supercell matrix &
              &must be integers.")
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
        tBadIntegratingKPoints = .true.
        call init(li1)
        call init(lr1)
        call getChildValue(value, "", 1, li1, 3, lr1)
        if (len(li1) < 1) then
          call detailedError(value, "At least one line must be specified.")
        end if
        allocate(tmpI1(len(li1)))
        allocate(kpts(3, 0:len(lr1)))
        call asVector(li1, tmpI1)
        call asArray(lr1, kpts(:,1:len(lr1)))
        kpts(:,0) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
        call destruct(li1)
        call destruct(lr1)
        if (any(tmpI1 < 0)) then
          call detailedError(value, "Interval steps must be greater equal to &
              &zero.")
        end if
        ctrl%nKPoint = sum(tmpI1)
        if (ctrl%nKPoint < 1) then
          call detailedError(value, "Sum of the interval steps must be greater &
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
          rTmp3 = (kpts(:,jj) - kpts(:,jj-1)) / real(tmpI1(jj), dp)
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
        tBadIntegratingKPoints = .false.

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
        call detailedError(value, "Invalid K-point scheme")
      end select
    end if

    if (ctrl%tSCC) then
      if (tBadIntegratingKPoints) then
        ii = 1
      else
        ii = 100
      end if
      call getChildValue(node, "MaxSCCIterations", ctrl%maxIter, ii)
    end if

    if (tBadIntegratingKPoints .and. ctrl%tSCC .and. ctrl%maxIter /= 1) then
      write(errorStr, "(A,I3)") "SCC cycle with these k-points probably will&
          & not correctly calculate many properties, SCC iterations set to:", ctrl%maxIter
      call warning(errorStr)
    end if
    if (tBadIntegratingKPoints .and. ctrl%tSCC .and. .not.ctrl%tReadChrg) then
      call warning("It is strongly suggested you use the ReadInitialCharges option.")
    end if

    call getChild(node, "OrbitalPotential", child, requested=.false.)
    if (.not. associated(child)) then
      ctrl%tDFTBU = .false.
      ctrl%DFTBUfunc = 0
    else
      call getChildValue(child, "Functional", buffer, "fll")
      select case(tolower(char(buffer)))
      case ("fll")
        ctrl%DFTBUfunc = plusUFunctionals%fll
      case ("psic")
        ctrl%DFTBUfunc = plusUFunctionals%pSic
      case default
        call detailedError(child,"Unknown orbital functional :"// char(buffer))
      end select

      allocate(ctrl%nUJ(geo%nSpecies))
      ctrl%nUJ = 0

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
        ctrl%nUJ(iSp1) = getLength(children)
        do ii = 1, ctrl%nUJ(iSp1)
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
          call convertByMul(char(modifier), energyUnits, child3, rTmp)
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
      end do

      do iSp1 = 1, geo%nSpecies
        ctrl%nUJ(iSp1) = len(lrN(iSp1))
      end do
      allocate(ctrl%UJ(maxval(ctrl%nUJ),geo%nSpecies))
      ctrl%UJ = 0.0_dp
      allocate(ctrl%niUJ(maxval(ctrl%nUJ),geo%nSpecies))
      ctrl%niUJ = 0
      do iSp1 = 1, geo%nSpecies
        call asArray(lrN(iSp1),ctrl%UJ(1:len(lrN(iSp1)),iSp1))
        allocate(iTmpN(len(liN(iSp1))))
        call asArray(liN(iSp1),iTmpN)
        ctrl%niUJ(1:len(liN(iSp1)),iSp1) = iTmpN(:)
        deallocate(iTmpN)
        call destruct(lrN(iSp1))
        call destruct(liN(iSp1))
      end do
      allocate(ctrl%iUJ(maxval(ctrl%niUJ),maxval(ctrl%nUJ),geo%nSpecies))
      ctrl%iUJ = 0
      do iSp1 = 1, geo%nSpecies
        do ii = 1, ctrl%nUJ(iSp1)
          allocate(iTmpN(ctrl%niUJ(ii,iSp1)))
          call get(li1N(iSp1),iTmpN,ii)
          ctrl%iUJ(1:ctrl%niUJ(ii,iSp1),ii,iSp1) = iTmpN(:)
          deallocate(iTmpN)
        end do
        call destruct(li1N(iSp1))
      end do

      deallocate(li1N)
      deallocate(lrN)
      deallocate(liN)

      ! sanity check time
      allocate(iTmpN(slako%orb%mShell))
      do iSp1 = 1, geo%nSpecies
        iTmpN = 0
        ! loop over number of blocks for that species
        do ii = 1, ctrl%nUJ(iSp1)
          iTmpN(ctrl%iUJ(1:ctrl%niUJ(ii,iSp1),ii,iSp1)) = &
              & iTmpN(ctrl%iUJ(1:ctrl%niUJ(ii,iSp1),ii,iSp1)) + 1
        end do
        if (any(iTmpN(:)>1)) then
          write(stdout, *)'Multiple copies of shells present in OrbitalPotential!'
          write(stdout, "(A,A3,A,I2)") &
              & 'The count for the occurance of shells of species ', &
              & trim(geo%speciesNames(iSp1)),' are:'
          write(stdout, *)iTmpN(1:slako%orb%nShell(iSp1))
          stop
        end if
      end do
      deallocate(iTmpN)

      ctrl%tDFTBU = .true.

    end if

    ! Dispersion
    call getChildValue(node, "Dispersion", value, "", child=child, &
        &allowEmptyValue=.true., dummyValue=.true.)
    if (associated(value)) then
      allocate(ctrl%dispInp)
      call readDispersion(child, geo, ctrl%dispInp)
    end if
    if (ctrl%tLatOpt .and. .not. geo%tPeriodic) then
      call error("Lattice optimization only applies for periodic structures.")
    end if

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
      if (ctrl%t3rd .and. ctrl%tOrbResolved) then
        call error("Only full third-order DFTB is compatible with orbital&
            & resolved SCC")
      end if
      if (ctrl%t3rd .or. ctrl%t3rdFull) then
        call getChild(node, 'HubbardDerivs', child, requested=.true.)
        allocate(ctrl%HubDerivs(slako%orb%mShell, geo%nSpecies))
        ctrl%hubDerivs(:,:) = 0.0_dp
        do iSp1 = 1, geo%nSpecies
          nShell = slako%orb%nShell(iSp1)
          if (ctrl%tOrbResolved) then
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
      end if
    end if

    call readDifferentiation(node, ctrl)

    if (ctrl%tSCC) then ! Force type
      call getChildValue(node, "ForceEvaluation", buffer, "Traditional", &
          & child=child)
      select case (tolower(unquote(char(buffer))))
      case("traditional")
        ctrl%forceType = 0
      case("dynamicst0")
        ctrl%forceType = 2
      case("dynamics")
        ctrl%forceType = 3
      case default
        call detailedError(child, "Invalid force evaluation method.")
      end select
    else
      ctrl%forceType = 0
    end if

    call readCustomisedHubbards(node, geo, slako%orb, ctrl%tOrbResolved, ctrl%hubbU)

  end subroutine readDFTBHam


  !> Reads inital charges
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
        call getChildValue(child2, "Atoms", buffer, child=child3, &
            &multiple=.true.)
        call convAtomRangeToInt(char(buffer), geo%speciesNames, &
            &geo%species, child3, pTmpI1)
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
        call getChildValue(child2, "Atoms", buffer, child=child3, &
            &multiple=.true.)
        call convAtomRangeToInt(char(buffer), geo%speciesNames, &
            &geo%species, child3, pTmpI1)
        call getChildValue(child2, "SpinPerAtom", rTmp)
        do jj = 1, size(pTmpI1)
          iAt = pTmpI1(jj)
          if (any(initSpins(:,iAt) /= 0.0_dp)) then
            call detailedWarning(child3, "Previoius setting for the spin &
                &of atom" // i2c(iAt) // " overwritten")
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
    type(control), intent(inout) :: ctrl


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
      ctrl%iDerivMethod = 1
      call getChildValue(val, "Delta", ctrl%deriv1stDelta, defDelta,&
          & modifier=modifier, child=child)
      call convertByMul(char(modifier), lengthUnits, child,&
          & ctrl%deriv1stDelta)
    case ("richardson")
      ctrl%iDerivMethod = 2
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
    type(control), intent(inout) :: ctrl

    type(fnode), pointer :: value, child, child2
    type(string) :: buffer
    real(dp) :: h5ScalingDef
    integer :: iSp

    ! X-H interaction corrections including H5 and damping
    ctrl%tDampH = .false.
    ctrl%h5SwitchedOn = .false.
    call getChildValue(node, "HCorrection", value, "None", child=child)
    call getNodeName(value, buffer)
    select case (char(buffer))
    case ("none")
      ! nothing to do
    case ("damping")
      ! Switch the correction on
      ctrl%tDampH = .true.
      call getChildValue(value, "Exponent", ctrl%dampExp)
    case ("h5")
      ! Switch the correction on
      ctrl%h5SwitchedOn = .true.

      call getChildValue(value, "RScaling", ctrl%h5RScale, 0.714_dp)
      call getChildValue(value, "WScaling", ctrl%h5WScale, 0.25_dp)

      allocate(ctrl%h5ElementPara(geo%nSpecies))
      call getChild(value, "H5Scaling", child2, requested=.false.)
      if (.not. associated(child2)) then
        call setChild(value, "H5scaling", child2)
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
        call getChildValue(child2, geo%speciesNames(iSp), ctrl%h5ElementPara(iSp), h5ScalingDef)
      end do
    case default
      call getNodeHSDName(value, buffer)
      call detailedError(child, "Invalid HCorrection '" // char(buffer) // "'")
    end select

  end subroutine readHCorrection


  !> Reads Slater-Koster files
  !> Should be replaced with a more sophisticated routine, once the new SK-format has been
  !> established
  subroutine readSKFiles(skFiles, nSpecies, slako, orb, angShells, orbRes, skInterMeth, repPoly)

    !> List of SK file names to read in for every interaction
    type(ListCharLc), intent(inout) :: skFiles(:,:)

    !> Nr. of species in the system
    integer, intent(in) :: nSpecies

    !> Data type for slako information
    type(slater), intent(inout) :: slako

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> For every species, a list of rank one arrays. Each array contains the angular momenta to pick
    !> from the appropriate SK-files.
    type(listIntR1), intent(inout) :: angShells(:)

    !> Are the Hubbard Us different for each l-shell?
    logical, intent(in) :: orbRes

    !> Method of the sk interpolation
    integer, intent(in) :: skInterMeth

    !> is this a polynomial or spline repulsive?
    logical, intent(in) :: repPoly(:,:)

    integer :: iSp1, iSp2, nSK1, nSK2, iSK1, iSK2, ind, nInt, iSh1
    integer :: angShell(maxL+1), nShell
    logical :: readRep, readAtomic
    character(lc) :: fileName
    real(dp), allocatable, target :: skHam(:,:), skOver(:,:)
    real(dp) :: dist
    type(TOldSKData), allocatable :: skData12(:,:), skData21(:,:)
    type(OSlakoEqGrid), allocatable :: pSlakoEqGrid1, pSlakoEqGrid2
    type(TRepSplineIn) :: repSplineIn1, repSplineIn2
    type(TRepPolyIn) :: repPolyIn1, repPolyIn2
    type(ORepSpline), allocatable :: pRepSpline
    type(ORepPoly), allocatable :: pRepPoly

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
    allocate(slako%repCont)
    call init(slako%repCont, nSpecies)

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
            write(stdout, "(2X,A)") trim(fileName)
            if (readRep .and. repPoly(iSp2, iSp1)) then
              call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, &
                  &repPolyIn=repPolyIn1)
            elseif (readRep) then
              call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic, &
                  &iSp1, iSp2, repSplineIn=repSplineIn1)
            else
              call readFromFile(skData12(iSK2,iSK1), fileName, readAtomic)
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
                    &repPolyIn=repPolyIn2)
              elseif (readRep) then
                call readFromFile(skData21(iSK1,iSK2), fileName, readAtomic, &
                    &iSp2, iSp1, repSplineIn=repSplineIn2)
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
        nInt = getNSKIntegrals(iSp1, iSp2, orb)
        allocate(skHam(size(skData12(1,1)%skHam, dim=1), nInt))
        allocate(skOver(size(skData12(1,1)%skOver, dim=1), nInt))
        call getFullTable(skHam, skOver, skData12, skData21, angShells(iSp1), &
            &angShells(iSp2))

        ! Add H/S tables to the containers for iSp1-iSp2
        dist = skData12(1,1)%dist
        allocate(pSlakoEqGrid1, pSlakoEqGrid2)
        call init(pSlakoEqGrid1, dist, skHam, skInterMeth)
        call init(pSlakoEqGrid2, dist, skOver, skInterMeth)
        call addTable(slako%skHamCont, pSlakoEqGrid1, iSp1, iSp2)
        call addTable(slako%skOverCont, pSlakoEqGrid2, iSp1, iSp2)
        deallocate(skHam)
        deallocate(skOver)
        if (iSp1 /= iSp2) then
          ! Heteronuclear interactions: the same for the reverse interaction
          allocate(skHam(size(skData12(1,1)%skHam, dim=1), nInt))
          allocate(skOver(size(skData12(1,1)%skOver, dim=1), nInt))
          call getFullTable(skHam, skOver, skData21, skData12, angShells(iSp2),&
              &angShells(iSp1))
          allocate(pSlakoEqGrid1, pSlakoEqGrid2)
          call init(pSlakoEqGrid1, dist, skHam, skInterMeth)
          call init(pSlakoEqGrid2, dist, skOver, skInterMeth)
          call addTable(slako%skHamCont, pSlakoEqGrid1, iSp2, iSp1)
          call addTable(slako%skOverCont, pSlakoEqGrid2, iSp2, iSp1)
          deallocate(skHam)
          deallocate(skOver)
        end if
        deallocate(skData12)
        deallocate(skData21)

        ! Add repulsives to the containers.
        if (repPoly(iSp2, iSp1)) then
          allocate(pRepPoly)
          call init(pRepPoly, repPolyIn1)
          call addRepulsive(slako%repCont, pRepPoly, iSp1, iSp2)
          deallocate(pRepPoly)
        else
          allocate(pRepSpline)
          call init(pRepSpline, repSplineIn1)
          call addRepulsive(slako%repCont, pRepSpline, iSp1, iSp2)
          deallocate(pRepSpline)
          deallocate(repSplineIn1%xStart)
          deallocate(repSplineIn1%spCoeffs)
        end if
        if (iSp1 /= iSp2) then
          if (repPoly(iSp1, iSp2)) then
            allocate(pRepPoly)
            call init(pRepPoly, repPolyIn2)
            call addRepulsive(slako%repCont, pRepPoly, iSp2, iSp1)
            deallocate(pRepPoly)
          else
            allocate(pRepSpline)
            call init(pRepSpline, repSplineIn2)
            call addRepulsive(slako%repCont, pRepSpline, iSp2, iSp1)
            deallocate(pRepSpline)
            deallocate(repSplineIn2%xStart)
            deallocate(repSplineIn2%spCoeffs)
          end if
        end if
      end do lpSp2
    end do lpSp1
    write(stdout, "(A)") "Done."

  end subroutine readSKFiles


  !> Checks if the provided set of SK-tables for a the interactions A-B and B-A are consistent.
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
    type(TRepSplineIn), intent(in) :: repIn1

    !> Repulsive spline for interaction B-A
    type(TRepSplineIn), intent(in) :: repIn2

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
    type(TRepPolyIn), intent(in) :: repIn1

    !> Repulsive polynomial for interaction B-A
    type(TRepPolyIn), intent(in) :: repIn2

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
  pure function getNSKIntegrals(sp1, sp2, orb) result(nInt)

    !> Index of the first species
    integer, intent(in) :: sp1

    !> Index of the second species
    integer, intent(in) :: sp2

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

    !> Nr. of Slater-Koster interactions
    integer :: nInt

    integer :: iSh1, iSh2

    nInt = 0
    do iSh1 = 1, orb%nShell(sp1)
      do iSh2 = 1, orb%nShell(sp2)
        nInt = nInt + min(orb%angShell(iSh2, sp2), orb%angShell(iSh1, sp1)) + 1
      end do
    end do

  end function getNSKIntegrals


  !> Creates from the columns of the Slater-Koster files for A-B and B-A a full table for A-B,
  !> containing all integrals.
  subroutine getFullTable(skHam, skOver, skData12, skData21, angShells1, &
      &angShells2)

    !> Resulting table of H integrals
    real(dp), intent(out) :: skHam(:,:)

    !> Resulting table of S integrals
    real(dp), intent(out) :: skOver(:,:)

    !> Contains all SK files describing interactions for A-B
    type(TOldSKData), intent(in), target :: skData12(:,:)

    !> Contains all SK files describing interactions for B-A
    type(TOldSKData), intent(in), target :: skData21(:,:)

    !> Angular momenta to pick from the SK-files for species A
    type(listIntR1), intent(inout) :: angShells1

    !> Angular momenta to pick from the SK-files for species B
    type(listIntR1), intent(inout) :: angShells2

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
    type(control), intent(inout) :: ctrl

    type(fnode), pointer :: child
    logical :: tWriteDetailedOutDef

  #:if WITH_SOCKETS
    tWriteDetailedOutDef = .not. allocated(ctrl%socketInput)
  #:else
    tWriteDetailedOutDef = .true.
  #:endif
    call getChildValue(node, "WriteDetailedOut", ctrl%tWriteDetailedOut, tWriteDetailedOutDef)

    call getChildValue(node, "WriteAutotestTag", ctrl%tWriteTagged, .false.)
    call getChildValue(node, "WriteDetailedXML", ctrl%tWriteDetailedXML, &
        &.false.)
    call getChildValue(node, "WriteResultsTag", ctrl%tWriteResultsTag, &
        &.false.)


    if (.not.(ctrl%tMD.or.ctrl%tGeoOpt)) then
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
    call getChildValue(node, "MinimiseMemoryUsage", ctrl%tMinMemory, .false., child=child)
    if (ctrl%tMinMemory) then
      call detailedWarning(child, "Memory minimisation is not working currently, normal calculation&
          & will be used instead")
    end if
    call getChildValue(node, "ShowFoldedCoords", ctrl%tShowFoldedCoord, .false.)
  #:if DEBUG > 0
    call getChildValue(node, "TimingVerbosity", ctrl%timingLevel, -1)
  #:else
    call getChildValue(node, "TimingVerbosity", ctrl%timingLevel, 0)
  #:endif

    if (ctrl%tReadChrg) then
      call getChildValue(node, "ReadChargesAsText", ctrl%tReadChrgAscii, .false.)
    end if
    call getChildValue(node, "WriteChargesAsText", ctrl%tWriteChrgAscii, .false.)

    ctrl%tSkipChrgChecksum = .false.
    if (.not. ctrl%tFixEf .and. ctrl%tReadChrg) then
      call getChildValue(node, "SkipChargeTest", ctrl%tSkipChrgChecksum, .false.)
    end if

  end subroutine readOptions


  !> Reads in dispersion related settings
  subroutine readDispersion(node, geo, input)

    !> Node to parse
    type(fnode), pointer :: node

    !> geometry, including atomic information
    type(TGeometry), intent(in) :: geo

    !> dispersion data on exit
    type(DispersionInp), intent(out) :: input

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
#:if WITH_DFTD3
      allocate(input%dftd3)
      call readDispDFTD3(dispModel, input%dftd3)
#:else
      call detailedError(node, "Program had been compiled without DFTD3 support")
#:endif
    case default
      call detailedError(node, "Invalid dispersion model name.")
    end select

  end subroutine readDispersion


  !> Reads in the dispersion input data for the Slater-Kirkwood dispersion modell.
  subroutine readDispSlaKirk(node, geo, input)

    !> Node to process
    type(fnode), pointer :: node

    !> Geometry of the current system
    type(TGeometry), intent(in) :: geo

    !> Contains the input for the dispersion module on exit
    type(DispSlaKirkInp), intent(out) :: input

    type(fnode), pointer :: value, value2, child, child2, child3
    type(string) :: buffer, modif, modif2, modifs(3)
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
    call getChildValue(node, "PolarRadiusCharge", value, child=child, &
        &modifier=modif)
    call getNodeName(value, buffer)
    select case (char(buffer))
    case (textNodeName)
      call getChildValue(child, "", tmpR2, modifier=modif)
      if (len(modif) > 0) then
        call splitModifier(char(modif), child, modifs)
        call convertByMul(char(modifs(1)), volumeUnits, child, tmpR2(1,:),&
            &.false.)
        call convertByMul(char(modifs(2)), lengthUnits, child, tmpR2(2,:),&
            &.false.)
        call convertByMul(char(modifs(3)), chargeUnits, child, tmpR2(3,:),&
            &.false.)
      end if

    case ("hybriddependentpol")
      if (len(modif) > 0) then
        call detailedError(child, "PolarRadiusCharge is not allowed to carry &
            &a modifier, if the HybridDependentPol method is used.")
      end if
      allocate(rCutoffs(geo%nSpecies))
      allocate(tmp2R2(13, geo%nSpecies))
      do iSp1 = 1, geo%nSpecies
        call getChildValue(value, geo%speciesNames(iSp1), value2, &
            &child=child2, dummyValue=.true.)
        call getChildValue(child2, "CovalentRadius", rCutoffs(iSp1), &
            &modifier=modif2, child=child3)
        call convertByMul(char(modif2), lengthUnits, child3, &
            &rCutoffs(iSp1))
        call getChildValue(child2, "HybridPolarisations", tmp2R2(:, iSp1), &
            &modifier=modif2, child=child3)
        if (len(modif2) > 0) then
          call splitModifier(char(modif2), child, modifs)
          call convertByMul(char(modifs(1)), volumeUnits, child, &
              &tmp2R2(1:6, iSp1), .false.)
          call convertByMul(char(modifs(2)), lengthUnits, child, &
              &tmp2R2(7:12, iSp1), .false.)
          call convertByMul(char(modifs(3)), chargeUnits, child, &
              &tmp2R2(13, iSp1), .false.)
        end if
      end do
      mCutoff = 2.0_dp * maxval(rCutoffs)
      if (geo%tPeriodic) then
        call getCellTranslations(cellVec, rCellVec, geo%latVecs, &
            & geo%recVecs2p, mCutoff)
      else
        allocate(cellVec(3, 1))
        allocate(rCellVec(3, 1))
        cellVec(:, 1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
        rCellVec(:, 1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
      end if
      call init(neighs, geo%nAtom, 10)
      if (geo%tPeriodic) then
        ! Make some guess for the nr. of all interacting atoms
        nAllAtom = int((real(geo%nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
      else
        nAllAtom = geo%nAtom
      end if
      allocate(coords(3, nAllAtom))
      allocate(img2CentCell(nAllAtom))
      allocate(iCellVec(nAllAtom))
      call updateNeighbourList(coords, img2CentCell, iCellVec, neighs, &
          &nAllAtom, geo%coords, mCutoff, rCellVec)
      allocate(nNeighs(geo%nAtom))
      nNeighs(:) = 0
      do iAt1 = 1, geo%nAtom
        iSp1 = geo%species(iAt1)
        do iNeigh = 1, neighs%nNeighbourSK(iAt1)
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
      call detailedError(value, "Invalid method for PolarRadiusCharge.")
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
    type(DispUffInp), intent(out) :: input

    type(string) :: buffer
    type(fnode), pointer :: child, value, child2
    integer :: iSp
    logical :: found

    call getChildValue(node, "Parameters", value, child=child)
    allocate(input%distances(geo%nSpecies))
    allocate(input%energies(geo%nSpecies))
    call getNodeName(value, buffer)
    select case(char(buffer))
    case("uffparameters")
      do iSp = 1, geo%nSpecies
        call getUffValues(geo%speciesNames(iSp), input%distances(iSp), &
            &input%energies(iSp), found)
        if (.not. found) then
          call detailedError(value, "UFF parameters for species '" // geo&
              &%speciesNames(iSp) // "' not found.")
        end if
      end do
    case default
      call setUnprocessed(value)
      do iSp = 1, geo%nSpecies
        call getChild(child, geo%speciesNames(iSp), child2)
        call getChildValue(child2, "Distance", input%distances(iSp), &
            &modifier=buffer)
        call convertByMul(char(buffer), lengthUnits, child, &
            &input%distances(iSp))
        call getChildValue(child2, "Energy", input%energies(iSp), &
            &modifier=buffer)
        call convertByMul(char(buffer), energyUnits, child, &
            &input%energies(iSp))
      end do
    end select

  end subroutine readDispVdWUFF

#:if WITH_DFTD3


  !> Reads in initialization data for the DFTD3 dispersion module.
  subroutine readDispDFTD3(node, input)

    !> Node to process.
    type(fnode), pointer :: node

    !> Filled input structure on exit.
    type(DispDftD3Inp), intent(out) :: input

    type(fnode), pointer :: child, childval
    type(string) :: buffer

    call getChildValue(node, "Damping", childval, default="BeckeJohnson", &
        & child=child)
    call getNodeName(childval, buffer)
    select case (char(buffer))
    case ("beckejohnson")
      input%tBeckeJohnson = .true.
      call getChildValue(childval, "a1", input%a1, default=0.5719_dp)
      call getChildValue(childval, "a2", input%a2, default=3.6017_dp)
      ! Alpha is not used in BJ-damping, however, there are unused terms,
      ! which are calculated with alpha nevertheless, so set the default
      ! as found in dftd3 code.
      input%alpha6 = 14.0_dp
    case ("zerodamping")
      input%tBeckeJohnson = .false.
      call getChildValue(childval, "sr6", input%sr6)
      ! Although according to the documentation, this parameter is not used
      ! when calculating zero damping, results do change, when this parameter
      ! is changed. We set it to the value found in dftd3 code
      input%sr8 = 1.0_dp
      call getChildValue(childval, "alpha6", input%alpha6, default=14.0_dp)
    case default
      call getNodeHSDName(childval, buffer)
      call detailedError(child, "Invalid damping method '" // char(buffer) // "'")
    end select
    call getChildValue(node, "s6", input%s6, default=1.0_dp)
    call getChildValue(node, "s8", input%s8, default=0.5883_dp)
    call getChildValue(node, "cutoff", input%cutoff, default=sqrt(9000.0_dp), &
        & modifier=buffer, child=child)
    call convertByMul(char(buffer), lengthUnits, child, input%cutoff)
    call getChildValue(node, "cutoffcn", input%cutoffCN, default=40.0_dp, &
        & modifier=buffer, child=child)
    call convertByMul(char(buffer), lengthUnits, child, input%cutoffCN)
    call getChildValue(node, "threebody", input%threebody, default=.false.)
    ! D3H5 - additional H-H repulsion
    call getChildValue(node, "hhrepulsion", input%hhrepulsion, default=.false.)

    input%numgrad = .false.

  end subroutine readDispDFTD3

#:endif


  !> reads in value of temperature for MD with sanity checking of the input
  subroutine readTemperature(node, ctrl)

    !> data to parse
    type(fnode), pointer :: node

    !> control data coming back
    type(control), intent(inout) :: ctrl

    type(string) :: modifier

    allocate(ctrl%tempSteps(1))
    allocate(ctrl%tempValues(1))
    allocate(ctrl%tempMethods(1))
    ctrl%tempMethods(1) = 1
    ctrl%tempSteps(1) = 1
    call getChildValue(node, "", ctrl%tempValues(1), modifier=modifier)
    call convertByMul(char(modifier), energyUnits, node, ctrl%tempValues(1))
    if (ctrl%tempValues(1) < 0.0_dp) then
      call detailedError(node, "Negative temperature.")
    end if
    if (ctrl%tempValues(1) < minTemp) then
      ctrl%tempValues(1) = minTemp
    end if

  end subroutine readTemperature


  !> reads a temperature profile for MD with sanity checking of the input
  subroutine readTemperatureProfile(node, modifier, ctrl)

    !> parser node contaning the relevant part of the user input
    type(fnode), pointer :: node

    !> unit modifier for the profile
    character(len=*), intent(in) :: modifier

    !> Control structure to populate
    type(control), intent(inout) :: ctrl


    !> Names of thermal profiles
    character(len=*), parameter :: tempMethodNames(3) = (/ 'constant   ', &
        &'linear     ', 'exponential' /)

    type(listString) :: ls
    type(listIntR1) :: li1
    type(listRealR1) :: lr1
    character(len=20), allocatable :: tmpC1(:)
    integer :: ii, jj

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
    lp2: do ii = 1, size(tmpC1)
      do jj = 1, size(tempMethodNames)
        if (trim(tmpC1(ii)) == tolower(trim(tempMethodNames(jj)))) then
          ctrl%tempMethods(ii) = jj
          cycle lp2
        end if
      end do
      call detailedError(node, "Invalid annealing method name '" &
          &// trim(tmpC1(ii)) // "'.")
    end do lp2

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

    call convertByMul(modifier, energyUnits, node, ctrl%tempValues)
    if (any(ctrl%tempValues < minTemp)) then
      ctrl%tempValues = max(ctrl%tempValues, minTemp)
    end if
    deallocate(tmpC1)

  end subroutine readTemperatureProfile


  !> Reads the excited state data block
  subroutine readExcited(node, ctrl)

    !> Node to parse
    type(fnode), pointer :: node

    !> Control structure to fill
    type(control), intent(inout) :: ctrl

    type(fnode), pointer :: child
  #:if WITH_ARPACK
    type(fnode), pointer :: child2
    type(string) :: buffer
    type(string) :: modifier
  #:endif

    ! Linear response stuff
    call getChild(node, "Casida", child, requested=.false.)

#:if not WITH_ARPACK

    if (associated(child)) then
      call detailedError(child, 'This DFTB+ binary has been compiled without support for linear&
          & response calculations (requires the ARPACK/ngARPACK library).')
    end if

#:else

    if (associated(child)) then

      ctrl%lrespini%tInit = .true.

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
      call convertByMul(char(modifier), energyUnits, child2, ctrl%lrespini%energyWindow)
      call getChildValue(child, "OscillatorWindow", ctrl%lrespini%oscillatorWindow, 0.0_dp, &
          & modifier=modifier,  child=child2)
      ctrl%lrespini%tOscillatorWindow = ctrl%lrespini%oscillatorWindow /= 0.0_dp
      call convertByMul(char(modifier), dipoleUnits, child2, ctrl%lrespini%oscillatorWindow)
      call getChildValue(child, "WriteMulliken", ctrl%lrespini%tMulliken, default=.false.)
      call getChildValue(child, "WriteCoefficients", ctrl%lrespini%tCoeffs, default=.false.)
      ctrl%lrespini%tGrndState = .false.
      if (ctrl%lrespini%tCoeffs) then
        call getChildValue(child, "TotalStateCoeffs", ctrl%lrespini%tGrndState, .false.)
      end if
      call getChildValue(child, "WriteEigenvectors", ctrl%lrespini%tPrintEigVecs, .false.)
      call getChildValue(child, "WriteXplusY", ctrl%lrespini%tXplusY, default=.false.)
      call getChildValue(child, "WriteSPTransitions", ctrl%lrespini%tSPTrans, default=.false.)
      call getChildValue(child, "WriteTransitions", ctrl%lrespini%tTrans, default=.false.)
      call getChildValue(child, "WriteTransitionDipole", ctrl%lrespini%tTradip, default=.false.)
      call getChildValue(child, "WriteStatusArnoldi", ctrl%lrespini%tArnoldi, default=.false.)
      call getChildValue(child, "TestArnoldi", ctrl%lrespini%tDiagnoseArnoldi, default=.false.)

    end if

#:endif

  end subroutine readExcited


  !> Reads the analysis block
  subroutine readAnalysis(node, ctrl, geo)

    !> Node to parse
    type(fnode), pointer :: node

    !> Control structure to fill
    type(control), intent(inout) :: ctrl

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    type(fnode), pointer :: val, child, child2, child3
    type(fnodeList), pointer :: children
    integer, allocatable :: pTmpI1(:)
    type(string) :: buffer
    integer :: nReg, iReg
    character(lc) :: strTmp
    type(listRealR1) :: lr1
    logical :: tPipekDense
    logical :: tWriteBandDatDef

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
        call getChildValue(child2, "Atoms", buffer, child=child3, &
            &multiple=.true.)
        call convAtomRangeToInt(char(buffer), geo%speciesNames, &
            &geo%species, child3, pTmpI1)
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

    call readElectrostaticPotential(node, geo, ctrl)

    call getChildValue(node, "MullikenAnalysis", ctrl%tPrintMulliken, .true.)
    call getChildValue(node, "AtomResolvedEnergies", ctrl%tAtomicEnergy, &
        &.false.)
    call getChildValue(node, "WriteEigenvectors", ctrl%tPrintEigVecs, .false.)

    if (ctrl%tPrintEigVecs .or. ctrl%lrespini%tPrintEigVecs) then
      call getChildValue(node, "EigenvectorsAsTxt", ctrl%tPrintEigVecsTxt, &
          & .false.)
    end if
  #:if WITH_SOCKETS
    tWriteBandDatDef = .not. allocated(ctrl%socketInput)
  #:else
    tWriteBandDatDef = .true.
  #:endif
    call getChildValue(node, "WriteBandOut", ctrl%tWriteBandDat, tWriteBandDatDef)
    call getChildValue(node, "CalculateForces", ctrl%tPrintForces, .false.)

  end subroutine readAnalysis


  !> Reads W values if required by settings in the Hamiltonian or the excited state
  subroutine readSpinConstants(hamNode, geo, slako, ctrl)

    !> node for Hamitonian data
    type(fnode), pointer :: hamNode

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Slater-Koster structure
    type(slater), intent(in) :: slako

    !> control structure
    type(control), intent(inout) :: ctrl

    type(fnode), pointer :: child
    logical :: tLRNeedsSpinConstants, tOrbResolvedW
    integer :: iSp1

    tLRNeedsSpinConstants = .false.

    if (ctrl%lrespini%tInit) then
      select case (ctrl%lrespini%sym)
      case ("T", "B", " ")
        tLRNeedsSpinConstants = .true.
      case ("S")
        tLRNeedsSpinConstants = .false.
      case default
      end select
    end if

    if (tLRNeedsSpinConstants .or. ctrl%tSpin) then
      allocate(ctrl%spinW(slako%orb%mShell, slako%orb%mShell, geo%nSpecies))
      ctrl%spinW(:,:,:) = 0.0_dp

      call getChild(hamNode, "SpinConstants", child)
      if (.not.ctrl%tOrbResolved) then
        call getChildValue(child, "ShellResolvedSpin", tOrbResolvedW, .false.)
      else
        tOrbResolvedW = .true.
      end if

      if (tOrbResolvedW) then
        ! potentially unique values for each shell
        do iSp1 = 1, geo%nSpecies
          call getChildValue(child, geo%speciesNames(iSp1),&
              & ctrl%spinW(:slako%orb%nShell(iSp1), :slako%orb%nShell(iSp1), iSp1))
        end do
      else
        ! only one value per atom
        do iSp1 = 1, geo%nSpecies
          call getChildValue(child, geo%speciesNames(iSp1),ctrl%spinW(1, 1, iSp1))
          ctrl%spinW(:slako%orb%nShell(iSp1), :slako%orb%nShell(iSp1), iSp1) =&
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


  !> Reads the parallel block.
  subroutine readParallel(root, parallelOpts)

    !> Root node eventually containing the current block
    type(fnode), pointer, intent(in) :: root

    !> Parallel settings
    type(TParallelOpts), allocatable, intent(out) :: parallelOpts

    type(fnode), pointer :: node

    call getChild(root, "Parallel", child=node, requested=.false.)
    if (withMpi .and. .not. associated(node)) then
      call setChild(root, "Parallel", node)
    end if
    if (associated(node)) then
      if (.not. withMpi) then
        call detailedWarning(node, "Settings will be read but ignored (compiled without MPI&
            & support)")
      end if
      allocate(parallelOpts)
      call getChildValue(node, "Groups", parallelOpts%nGroup, 1)
      call getChildValue(node, "UseOmpThreads", parallelOpts%tOmpThreads, .not. withMpi)
      call readBlacs(node, parallelOpts%blacsOpts)
    end if

  end subroutine readParallel


  !> Reads the blacs block.
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


  subroutine readElectrostaticPotential(node, geo, ctrl)

    !> Node containing optional electrostatic settings
    type(fnode), pointer, intent(in) :: node

    !> geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Control structure
    type(control), intent(inout) :: ctrl

    type(fnode), pointer :: child, child2, child3
    type(string) :: buffer, modifier
    type(listRealR1) :: lr1

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
    if (ctrl%tGeoOpt .or. ctrl%tMD) then
      call getChildValue(child, "AppendFile", ctrl%elStatPotentialsInp%tAppendEsp, .false.)
    end if
    call init(lr1)
    ! discrete points
    call getChildValue(child, "Points", child2, "", child=child3, &
        & modifier=modifier, allowEmptyValue=.true.)
    call getNodeName2(child2, buffer)
    if (char(buffer) /= "") then
      call getChildValue(child3, "", 3, lr1, modifier=modifier)
      allocate(ctrl%elStatPotentialsInp%espGrid(3,len(lr1)))
      call asArray(lr1, ctrl%elStatPotentialsInp%espGrid)
      if (geo%tPeriodic .and. (char(modifier) == "F" .or. char(modifier) == "f")) then
        ctrl%elStatPotentialsInp%espGrid = matmul(geo%latVecs, ctrl%elStatPotentialsInp%espGrid)
      else
        call convertByMul(char(modifier), lengthUnits, child3,&
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
    call convertByMul(char(modifier), lengthUnits, child2, ctrl%elStatPotentialsInp%softenEsp)

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
      call convertByMul(char(modifier), lengthUnits, child, r3Tmp)
      call convertByMul(char(modifier), lengthUnits, child, r3Tmpb)
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

end module parser
