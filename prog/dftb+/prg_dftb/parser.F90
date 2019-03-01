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
  use inputconversion
  use oldcompat
  use lapackroutines, only : matinv
  use periodic
  use simplealgebra, only: cross3, determinant33
  use dispersions
  use dftbplusu
  use slakocont
  use slakoeqgrid
  use repcont
  use repspline
  use reppoly
  use commontypes
  use oldskdata
  use timeprop_module
  use xmlf90
  use dftbp_forcetypes, only : forceTypes
  use mixer, only : mixerTypes
  use geoopt, only : geoOptTypes
#:if WITH_SOCKETS
  use ipisocket, only : IPI_PROTOCOLS
#:endif
  use elsiiface
  use elecsolvers, only : electronicSolverTypes
  use wrappedintrinsics
#:if WITH_TRANSPORT
  use poisson_init
  use libnegf_vars
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

    call readHamiltonian(hamNode, input%ctrl, input%geom, input%slako, input%transpar,&
        & input%ginfo%greendens, input%poisson)

    call getChild(root, "Dephasing", child, requested=.false.)
    if (associated(child)) then
      call detailedError(child, "Be patient... Dephasing feature will be available soon!")
      !call readDephasing(child, input%slako%orb, input%geom, input%transpar, input%ginfo%tundos)
    end if

  #:else

    if (associated(child)) then
      call detailedError(child, "Program had been compiled without transport enabled")
    end if

    call readHamiltonian(hamNode, input%ctrl, input%geom, input%slako)

  #:endif

    ! Geometry driver
    call getChildValue(root, "Driver", tmp, "", child=child, allowEmptyValue=.true.)
  #:if WITH_TRANSPORT
    call readDriver(tmp, child, input%geom, input%ctrl, input%transpar)
  #:else
    call readDriver(tmp, child, input%geom, input%ctrl)
  #:endif

    ! Analysis of properties
    call getChildValue(root, "Analysis", dummy, "", child=child, list=.true., &
        & allowEmptyValue=.true., dummyValue=.true.)

  #:if WITH_TRANSPORT
    call readAnalysis(child, input%ctrl, input%geom, input%slako%orb, input%transpar, &
        & input%ginfo%tundos)

    call finalizeNegf(input)
  #:else
    call readAnalysis(child, input%ctrl, input%geom, input%slako%orb)
  #:endif

    ! excited state options
    call getChildValue(root, "ExcitedState", dummy, "", child=child, list=.true., &
        & allowEmptyValue=.true., dummyValue=.true.)
    call readExcited(child, input%ctrl)

    ! Options for calculation
    call getChildValue(root, "Options", dummy, "", child=child, list=.true., &
        & allowEmptyValue=.true., dummyValue=.true.)
    call readOptions(child, input%ctrl)

    ! Read W values if needed by Hamitonian or excited state calculation
    call readSpinConstants(hamNode, input%geom, input%slako, input%ctrl)

    call readParallel(root, input)

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

    type(fnode), pointer :: value1, child
    type(string) :: buffer

    call getChildValue(node, "", value1, child=child)
    call getNodeName(value1, buffer)
    select case (char(buffer))
    case ("genformat")
      call readTGeometryGen(value1, input%geom)
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

    !> Control structure to be filled
    type(TGeometry), intent(in) :: geom

    !> Nr. of atoms in the system
    type(control), intent(inout) :: ctrl

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(in) :: transpar
  #:endif

    type(fnode), pointer :: child, child2, child3, value1, value2, field

    type(string) :: buffer, buffer2, modifier
  #:if WITH_SOCKETS
    character(lc) :: sTmp
  #:endif

    ! range of default atoms to move
    character(mc) :: atomsRange

    atomsRange = "1:-1"
  #:if WITH_TRANSPORT
    if (transpar%defined) then
      ! only those atoms in the device region
      write(atomsRange,"(I0,':',I0)")transpar%idxdevice
    end if
  #:endif

    ctrl%tGeoOpt = .false.
    ctrl%tCoordOpt = .false.
    ctrl%tLatOpt = .false.

    ctrl%iGeoOpt = geoOptTypes%none
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

      ctrl%iGeoOpt = geoOptTypes%steepestDesc
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
      call getChildValue(node, "MovedAtoms", buffer2, trim(atomsRange), child=child, &
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

      ctrl%iGeoOpt = geoOptTypes%conjugateGrad
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
      call getChildValue(node, "MovedAtoms", buffer2, trim(atomsRange), child=child, &
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

      ctrl%iGeoOpt = geoOptTypes%diis
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
      call getChildValue(node, "MovedAtoms", buffer2, trim(atomsRange), child=child, &
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

      ctrl%iGeoOpt = geoOptTypes%lbfgs

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
      call getChildValue(node, "MovedAtoms", buffer2, trim(atomsRange), child=child, &
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
      call getChildValue(node, "Atoms", buffer2, trim(atomsRange), child=child, &
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
      call getChildValue(node, "MovedAtoms", buffer2, trim(atomsRange), child=child, &
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

      call getChildValue(node, "Thermostat", value1, child=child)
      call getNodeName(value1, buffer2)

      call getChildValue(node, "ConvergentForcesOnly", ctrl%tConvrgForces, &
          & .true.)

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
            call convertByMul(char(modifier), timeUnits, child3, &
                & ctrl%wvScale)
            ctrl%wvScale = ctrl%deltaT / ctrl%wvScale
          else
            call error("Either CouplingStrength or Timescale must be set&
                & for Berendsen thermostats.")
          end if
        end if

        call getChildValue(value1, "AdaptFillingTemp", ctrl%tSetFillingTemp, &
            &.false.)

      case ("nosehoover")
        ctrl%iThermostat = 3
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

        call getChildValue(value1, "CouplingStrength", ctrl%wvScale, &
            & modifier=modifier, child=field)
        call convertByMul(char(modifier), freqUnits, field, ctrl%wvScale)

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

        call getChildValue(value1, "AdaptFillingTemp", ctrl%tSetFillingTemp, &
            &.false.)

      case ("andersen")
        ctrl%iThermostat = 1
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

        call getChildValue(value1, "ReselectProbability", ctrl%wvScale, &
            &child=child3)
        if (ctrl%wvScale <= 0.0_dp .or. ctrl%wvScale > 1.0_dp) then
          call detailedError(child3, &
              &"ReselectProbability must be in the range (0,1]!")
        end if
        call getChildValue(value1, "ReselectIndividually", ctrl%tRescale)
        call getChildValue(value1, "AdaptFillingTemp", ctrl%tSetFillingTemp, &
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
          call getChildValue(value1, "InitialTemperature", ctrl%tempAtom, &
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
        call getNodeHSDName(value1, buffer2)
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
      if (input%maxSCCIter <= 0) then
        call detailedError(pRoot,"MaxSccIterations must be >= 1");
      end if
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

    type(fnode), pointer :: value1, child
    type(string) :: buffer
    type(listIntR1) :: intBuffer
    type(listRealR1) :: realBuffer

    call getChildValue(node, "Constraints", value1, "", child=child, &
        &allowEmptyValue=.true.)
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
    type(control), intent(inout) :: ctrl

    !> Total number of all atoms
    integer, intent(in) :: nAtom

    type(fnode), pointer :: value1, child
    type(string) :: buffer, modifier
    type(listRealR1) :: realBuffer
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
#:if WITH_TRANSPORT
  subroutine readHamiltonian(node, ctrl, geo, slako, tp, greendens, poisson)
#:else
  subroutine readHamiltonian(node, ctrl, geo, slako)
#:endif

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(control), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Slater-Koster structure to be filled
    type(slater), intent(inout) :: slako

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson
  #:endif

    type(string) :: buffer

    call getNodeName(node, buffer)
    select case (char(buffer))
    case ("dftb")
  #:if WITH_TRANSPORT
      call readDFTBHam(node, ctrl, geo, slako, tp, greendens, poisson)
  #:else
      call readDFTBHam(node, ctrl, geo, slako)
  #:endif
    case default
      call detailedError(node, "Invalid Hamiltonian")
    end select

  end subroutine readHamiltonian

  !> Reads DFTB-Hamiltonian
#:if WITH_TRANSPORT
  subroutine readDFTBHam(node, ctrl, geo, slako, tp, greendens, poisson)
#:else
  subroutine readDFTBHam(node, ctrl, geo, slako)
#:endif

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(control), intent(inout) :: ctrl

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Slater-Koster structure to be filled
    type(slater), intent(inout) :: slako

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout)  :: tp

    !> Green's function paramenters
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

    !> Poisson solver paramenters
    type(TPoissonInfo), intent(inout) :: poisson
  #:endif

    type(fnode), pointer :: value1, value2, child, child2, child3, field
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
    real(dp) :: rSKCutOff

    ! Read in maximal angular momenta or selected shells
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
            do kk = 1, size(orbitalNames)
              if (tmpCh == trim(orbitalNames(kk))) then
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
        call getNodeHSDName(value1, buffer)
        call detailedError(child2, "Invalid shell specification method '" //&
            & char(buffer) // "'")
      end select
    end do

    ! Orbitals and angular momenta for the given shells (once the SK files contain the full
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
        call detailedError(value1, "Assymetric definition (both A-B and B-A must&
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

    call getChild(node, "TruncateSKRange", child, requested=.false.)
    if (associated(child)) then
      call warning("Artificially truncating the SK table, this is normally a bad idea!")
      call SKTruncations(child, rSKCutOff, skInterMeth)
      call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tOrbResolved,&
          & skInterMeth, repPoly, rSKCutOff)
    else
      rSKCutOff = 0.0_dp
      call readSKFiles(skFiles, geo%nSpecies, slako, slako%orb, angShells, ctrl%tOrbResolved,&
          & skInterMeth, repPoly)
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
    call getChildValue(node, "SCC", ctrl%tSCC, .false.)
    ifSCC: if (ctrl%tSCC) then
      call getChildValue(node, "ReadInitialCharges", ctrl%tReadChrg, .false.)
      if (.not. ctrl%tReadChrg) then
        call getInitialCharges(node, geo, ctrl%initialCharges)
      end if
      call getChildValue(node, "SCCTolerance", ctrl%sccTol, 1.0e-5_dp)

      ! temporararily removed until debugged
      !call getChildValue(node, "WriteShifts", ctrl%tWriteShifts, .false.)
      ctrl%tWriteShifts = .false.

      call getChildValue(node, "Mixer", value1, "Broyden", child=child)
      call getNodeName(value1, buffer)
      select case(char(buffer))

      case ("broyden")

        ctrl%iMixSwitch = mixerTypes%broyden
        call getChildValue(value1, "MixingParameter", ctrl%almix, 0.2_dp)
        call getChildValue(value1, "InverseJacobiWeight", ctrl%broydenOmega0, &
            &0.01_dp)
        call getChildValue(value1, "MinimalWeight", ctrl%broydenMinWeight, &
            &1.0_dp)
        call getChildValue(value1, "MaximalWeight", ctrl%broydenMaxWeight, &
            &1.0e5_dp)
        call getChildValue(value1, "WeightFactor", ctrl%broydenWeightFac, &
            &1.0e-2_dp)

      case ("anderson")
        ctrl%iMixSwitch = mixerTypes%anderson
        call getChildValue(value1, "MixingParameter", ctrl%almix, 0.05_dp)
        call getChildValue(value1, "Generations", ctrl%iGenerations, 4)
        call getChildValue(value1, "InitMixingParameter", ctrl%andersonInitMixing, 0.01_dp)
        call getChildValue(value1, "DynMixingParameters", value2, "", &
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
        call getChildValue(value1, "DiagonalRescaling", ctrl%andersonOmega0, &
            &1.0e-2_dp)

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

      if (geo%tPeriodic) then
        call getChildValue(node, "EwaldParameter", ctrl%ewaldAlpha, 0.0_dp)
        call getChildValue(node, "EwaldTolerance", ctrl%tolEwald, 1.0e-9_dp)
      end if

      ctrl%tMulliken = .true.
      call readHCorrection(node, geo, ctrl)

    end if ifSCC

    ! Spin calculation
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

#:if WITH_TRANSPORT
    if (ctrl%tSpin .and. tp%ncont > 0) then
       call detailedError(child, "Spin-polarized transport is under development" //&
             & "and not currently available")
    end if
#:endif

    ! temporararily removed until debugged
    !if (.not. ctrl%tscc) then
    !  !! In a non-SCC calculation it is possible to upload charge shifts
    !  !! This is useful if the calculation can jump directly to the Analysis block
    !  call getChildValue(node, "ReadShifts", ctrl%tReadShifts, .false.)
    !end if
    ctrl%tReadShifts = .false.

    ! External electric field
    call getChildValue(node, "ElectricField", value1, "", child=child, &
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
        call getChildValue(child2, "CoordsAndCharges", value1, &
            &modifier=modifier, child=child3)
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
          open(newunit=fp, file=unquote(char(buffer2)), form="formatted", status="old",&
              & action="read", iostat=iErr)
          if (iErr /= 0) then
            call detailedError(value1, "Could not open file '" &
                &// trim(unquote(char(buffer2))) // "' for direct reading" )
          end if
          read(fp, *, iostat=iErr) tmpR2
          if (iErr /= 0) then
            call detailedError(value1, "Error during direct reading '" &
                &// trim(unquote(char(buffer2))) // "'")
          end if
          close(fp)
          ctrl%nExtChrg = ctrl%nExtChrg + ind
        case default
          call detailedError(value1, "Invalid block name")
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

    ! Filling (temperature only read, if AdaptFillingTemp was not set for the selected MD
    ! thermostat.)
    call getChildValue(node, "Filling", value1, "Fermi", child=child)
    call getNodeName(value1, buffer)

    select case (char(buffer))
    case ("fermi")
      ctrl%iDistribFn = 0 ! Fermi function
    case ("methfesselpaxton")
      ! Set the order of the Methfessel-Paxton step function approximation, defaulting to 2
      call getChildValue(value1, "Order", ctrl%iDistribFn, 2)
      if (ctrl%iDistribFn < 1) then
        call getNodeHSDName(value1, buffer)
        write(errorStr, "(A,A,A,I4)")"Unsuported filling mode for '", &
            & char(buffer),"' :",ctrl%iDistribFn
        call detailedError(child, errorStr)
      end if
    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Invalid filling method '" //char(buffer)// "'")
    end select

    if (.not. ctrl%tSetFillingTemp) then
      call getChildValue(value1, "Temperature", ctrl%tempElec, 0.0_dp, &
          &modifier=modifier, child=field)
      call convertByMul(char(modifier), energyUnits, field, ctrl%tempElec)
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
      call convertByMul(char(modifier), energyUnits, child3, ctrl%Ef)
    end if

    if (geo%tPeriodic .and. .not.ctrl%tFixEf) then
      call getChildValue(value1, "IndependentKFilling", ctrl%tFillKSep, .false.)
    end if

    ! Electronic solver
    call getChildValue(node, "Eigensolver", value1, "RelativelyRobust")
    call getNodeName(value1, buffer)

    select case(char(buffer))

    case ("qr")
      ctrl%solver%isolver = electronicSolverTypes%qr

    case ("divideandconquer")
      ctrl%solver%isolver = electronicSolverTypes%divideandconquer

    case ("relativelyrobust")
      ctrl%solver%isolver = electronicSolverTypes%relativelyrobust

    case ("elpa")
      ctrl%solver%isolver = electronicSolverTypes%elpa
      allocate(ctrl%solver%elsi)
      ctrl%solver%elsi%iSolver = ctrl%solver%isolver
      call getChildValue(value1, "Mode", ctrl%solver%elsi%elpaSolver, 2)

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
      call getChildValue(value1, "Poles", ctrl%solver%elsi%pexsiNPole, 20)
      call getChildValue(value1, "ProcsPerPole", ctrl%solver%elsi%pexsiNpPerPole, 1)
      call getChildValue(value1, "muPoints", ctrl%solver%elsi%pexsiNMu, 2)
      call getChildValue(value1, "SymbolicFactorProcs", ctrl%solver%elsi%pexsiNpSymbo, 1)
      call getChildValue(value1, "SpectralRadius", ctrl%solver%elsi%pexsiDeltaE, 10.0_dp,&
          & modifier=modifier, child=child)
      call convertByMul(char(modifier), energyUnits, child, ctrl%solver%elsi%pexsiDeltaE)

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
      if (tp%defined .and. .not.tp%taskUpload) then
        call detailederror(node, "greensfunction solver cannot be used "// &
            &  "when task = contactHamiltonian")
      end if
      call readGreensFunction(value1, greendens, tp, ctrl%tempElec)
      ! fixEf also avoids checks of total charge in initQFromFile
      ctrl%tFixEf = .true.
      if (geo%tPeriodic .and. greendens%doLocalCurr) then
         call detailedError(value1, "Local Currents in periodic systems still needs" //&
              " debugging and will be available soon")
      end if
    case ("transportonly")
      if (ctrl%tGeoOpt .or. ctrl%tMD) then
        call detailederror(node, "transportonly cannot be used with relaxations or md")
      end if
      if (tp%defined .and. .not.tp%taskUpload) then
        call detailederror(node, "transportonly cannot be used when "// &
            &  "task = contactHamiltonian")
      end if
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
      call getChildValue(value1, "Sparse", ctrl%solver%elsi%elsiCsr, .false.)
      if (ctrl%t2Component) then
        call detailedError(value1,"Two-component hamiltonians currently cannot be used with sparse&
            & ELSI solvers")
      end if
    end if

  #:if WITH_TRANSPORT
    if (all(ctrl%solver%isolver /= [electronicSolverTypes%GF,electronicSolverTypes%OnlyTransport])&
        & .and. tp%taskUpload) then
      call detailedError(value1, "Eigensolver incompatible with transport calculation&
          & (GreensFunction or TransportOnly required)")
    end if
  #:endif

    ! Charge
    call getChildValue(node, "Charge", ctrl%nrChrg, 0.0_dp)

    ! Assume SCC can has usual default number of steps if needed
    tBadIntegratingKPoints = .false.

    ! K-Points
    if (geo%tPeriodic) then
      call getChildValue(node, "KPointsAndWeights", value1, child=child, &
          &modifier=modifier)
      call getNodeName(value1, buffer)
      select case(char(buffer))

      case ("supercellfolding")
        tBadIntegratingKPoints = .false.
        if (len(modifier) > 0) then
          call detailedError(child, "No modifier is allowed, if the &
              &SupercellFolding scheme is used.")
        end if
        call getChildValue(value1, "", coeffsAndShifts)
        if (abs(determinant33(coeffsAndShifts(:,1:3))) - 1.0_dp < -1e-6_dp) then
          call detailedError(value1, "Determinant of the supercell matrix must &
              &be greater than 1")
        end if
        if (any(abs(modulo(coeffsAndShifts(:,1:3) + 0.5_dp, 1.0_dp) - 0.5_dp) &
            &> 1e-6_dp)) then
          call detailedError(value1, "The components of the supercell matrix &
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
        call detailedError(value1, "Invalid K-point scheme")
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
      call readDispersion(child, geo, ctrl%dispInp)
    end if
    if (ctrl%tLatOpt .and. .not. geo%tPeriodic) then
      call error("Lattice optimization only applies for periodic structures.")
    end if

    ctrl%tPoisson = .false.

  #:if WITH_TRANSPORT
    ! Read in which kind of electrostatics method to use.
    call getChildValue(node, "Electrostatics", value1, "GammaFunctional", &
        &child=child)
    call getNodeName(value1, buffer)
    select case (char(buffer))
    case ("gammafunctional")
      if (tp%taskUpload .and. ctrl%tSCC) then
        call detailedError(child, "GammaFunctional not available, if you upload contacts in an SCC&
            & calculation.")
      end if
    case ("poisson")
      ctrl%tPoisson = .true.
      call readPoisson(value1, poisson, geo%tPeriodic, tp%tPeriodic1D)
    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Unknown electrostatics '" // char(buffer) // "'")
    end select
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
        ctrl%forceType = forceTypes%orig
      case("dynamicst0")
        ctrl%forceType = forceTypes%dynamicT0
      case("dynamics")
        ctrl%forceType = forceTypes%dynamicTFinite
      case default
        call detailedError(child, "Invalid force evaluation method.")
      end select
    else
      ctrl%forceType = forceTypes%orig
    end if

    call readCustomisedHubbards(node, geo, slako%orb, ctrl%tOrbResolved, ctrl%hubbU)

  end subroutine readDFTBHam


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
    call convertByMul(char(modifier), lengthUnits, field, truncationCutOff)

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

    type(fnode), pointer :: value1, child, child2
    type(string) :: buffer
    real(dp) :: h5ScalingDef
    integer :: iSp

    ! X-H interaction corrections including H5 and damping
    ctrl%tDampH = .false.
    ctrl%h5SwitchedOn = .false.
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
      ! Switch the correction on
      ctrl%h5SwitchedOn = .true.

      call getChildValue(value1, "RScaling", ctrl%h5RScale, 0.714_dp)
      call getChildValue(value1, "WScaling", ctrl%h5WScale, 0.25_dp)

      allocate(ctrl%h5ElementPara(geo%nSpecies))
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
        call getChildValue(child2, geo%speciesNames(iSp), ctrl%h5ElementPara(iSp), h5ScalingDef)
      end do
    case default
      call getNodeHSDName(value1, buffer)
      call detailedError(child, "Invalid HCorrection '" // char(buffer) // "'")
    end select

  end subroutine readHCorrection


  !> Reads Slater-Koster files
  !> Should be replaced with a more sophisticated routine, once the new SK-format has been
  !> established
  subroutine readSKFiles(skFiles, nSpecies, slako, orb, angShells, orbRes, skInterMeth, repPoly,&
      & truncationCutOff)

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

    !> Distances to artificially truncate tables of SK integrals
    real(dp), intent(in), optional :: truncationCutOff

    integer :: iSp1, iSp2, nSK1, nSK2, iSK1, iSK2, ind, nInteract, iSh1
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
        nInteract = getNSKIntegrals(iSp1, iSp2, orb)
        allocate(skHam(size(skData12(1,1)%skHam, dim=1), nInteract))
        allocate(skOver(size(skData12(1,1)%skOver, dim=1), nInteract))
        call getFullTable(skHam, skOver, skData12, skData21, angShells(iSp1), &
            &angShells(iSp2))

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
              &angShells(iSp1))
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
    call getChildValue(node, "TimingVerbosity", ctrl%timingLevel, 1)
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

    type(fnode), pointer :: value1, value2, child, child2, child3
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
    call getChildValue(node, "PolarRadiusCharge", value1, child=child, &
        &modifier=modif)
    call getNodeName(value1, buffer)
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
        call getChildValue(value1, geo%speciesNames(iSp1), value2, &
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
    type(DispUffInp), intent(out) :: input

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
      call getChildValue(child, "WriteSPTransitions", ctrl%lrespini%tSPTrans, default=.false.)
      call getChildValue(child, "WriteTransitions", ctrl%lrespini%tTrans, default=.false.)
      call getChildValue(child, "WriteTransitionDipole", ctrl%lrespini%tTradip, default=.false.)
      call getChildValue(child, "WriteStatusArnoldi", ctrl%lrespini%tArnoldi, default=.false.)
      call getChildValue(child, "TestArnoldi", ctrl%lrespini%tDiagnoseArnoldi, default=.false.)

      if (ctrl%tForces .or. ctrl%tPrintForces) then
        call getChildValue(child, "ExcitedStateForces", ctrl%tCasidaForces, default=.true.)
      end if

    end if

  #:endif

  end subroutine readExcited


  !> Reads the analysis block
#:if WITH_TRANSPORT
  subroutine readAnalysis(node, ctrl, geo, orb, transpar, tundos)
#:else
  subroutine readAnalysis(node, ctrl, geo, orb)
#:endif

    !> Node to parse
    type(fnode), pointer :: node

    !> Control structure to fill
    type(control), intent(inout) :: ctrl

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Orbital
    type(TOrbitals), intent(in) :: orb

  #:if WITH_TRANSPORT
    !> Transport parameters
    type(TTransPar), intent(inout) :: transpar

    !> Tunneling and Dos parameters
    type(TNEGFTunDos), intent(inout) :: tundos
  #:endif

    type(fnode), pointer :: val, child, child2, child3
    type(fnodeList), pointer :: children
    integer, allocatable :: pTmpI1(:)
    type(string) :: buffer
    integer :: nReg, iReg
    character(lc) :: strTmp
    type(listRealR1) :: lr1
    logical :: tPipekDense
    logical :: tWriteBandDatDef, tHaveEigenDecomposition

    tHaveEigenDecomposition = .false.
    if (any(ctrl%solver%isolver == [electronicSolverTypes%qr,&
        & electronicSolverTypes%divideandconquer, electronicSolverTypes%relativelyrobust,&
        & electronicSolverTypes%elpa])) then
      tHaveEigenDecomposition = .true.
    end if

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

    end if

    ! Is this compatible with Poisson solver use?
    call readElectrostaticPotential(node, geo, ctrl)

    call getChildValue(node, "MullikenAnalysis", ctrl%tPrintMulliken, .true.)
    call getChildValue(node, "AtomResolvedEnergies", ctrl%tAtomicEnergy, &
        &.false.)

    call getChildValue(node, "CalculateForces", ctrl%tPrintForces, .false.)

    call getChild(node, "ElectronDynamics", child=child, requested=.false.)
    if (associated(child)) then
       allocate(ctrl%elecDynInp)
       call readElecDynamics(child, ctrl%elecDynInp, geo, ctrl%masses)
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
      call readTunAndDos(child, orb, geo, tundos, transpar, ctrl%tempElec)
    endif
  #:endif


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


  !> Reads the electron dynamics block
  subroutine readElecDynamics(node, input, geo, masses)

    !> input data to parse
    type(fnode), pointer :: node
    type(TGeometry), intent(in) :: geo
    !> masses to be returned
    real(dp), allocatable, intent(inout) :: masses(:)

    !> ElecDynamicsInp instance
    type(TElecDynamicsInp), intent(inout) :: input

    type(fnode), pointer :: value1, value2, child, child2
    type(string) :: buffer, buffer2, modifier
    logical :: ppRangeInvalid
    real (dp) :: defPpRange(2)

  #:if WITH_MPI
    if (associated(node)) then
       call detailedError(node, 'This DFTB+ binary has been compiled with MPI settings and &
            & electron dynamics are currently not supported.')
    end if
  #:endif

    call getChildValue(node, "Steps", input%steps)
    call getChildValue(node, "TimeStep", input%dt, modifier=modifier, &
         & child=child)
    call convertByMul(char(modifier), timeUnits, child, input%dt)

    call getChildValue(node, "FieldStrength", input%tdfield, modifier=modifier, child=child)
    call convertByMul(char(modifier), EFieldUnits, child, input%tdfield)

    call getChildValue(node, "Populations", input%tPopulations, .false.)
    call getChildValue(node, "WriteFrequency", input%writeFreq, 50)
    call getChildValue(node, "Restart", input%tRestart, .false.)
    call getChildValue(node, "WriteRestart", input%tWriteRestart, .true.)
    call getChildValue(node, "RestartFrequency", input%restartFreq, input%Steps / 10)
    call getChildValue(node, "Forces", input%tForces, .false.)
!    call getChildValue(node, "WritePairWiseEnergy", input%tPairWise, .false.)
    call getChildValue(node, "WriteBondEnergy", input%tBondE, .false.)
    call getChildValue(node, "WriteBondOrder", input%tBondO, .false.)
    call getChildValue(node, "OnsiteGradients", input%tOnsiteGradients, .false.)
    call getChildValue(node, "Pump", input%tPump, .false.)

    if (input%tPump) then
      call getChildValue(node, "PumpProbeFrames", input%tdPPFrames)
      defPpRange = [0.0_dp, input%steps * input%dt]
      call getChildValue(node, "PumpProbeRange", input%tdPpRange, defPprange, modifier=modifier,&
           & child=child)
       call convertByMul(char(modifier), timeUnits, child, input%tdPpRange)

      ppRangeInvalid = (input%tdPpRange(2) <= input%tdPpRange(1))&
           & .or. (input%tdPprange(1) < defPpRange(1)) .or. (input%tdPpRange(2) > defPpRange(2))
      if (ppRangeInvalid) then
         call detailederror(child, "Wrong definition of PumpProbeRange")
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

    !! Different perturbation types
    call getChildValue(node, "Perturbation", value1, "None", child=child)
    call getNodeName(value1, buffer)
    select case(char(buffer))

    case ("kick")
       input%pertType = iKick
       call getChildValue(value1, "PolarizationDirection", input%polDir)
       if (input%polDir < 1 .or. input%polDir > 4) then
          call detailedError(child, "Wrong specified polarization direction")
       end if
       call getChildValue(value1, "SpinType", buffer2, "Singlet")
       select case(unquote(char(buffer2)))
       case ("singlet", "Singlet")
          input%spType = iTDSinglet
       case ("triplet", "Triplet")
          input%spType = iTDTriplet
       case default
          call detailedError(value1, "Unknown spectrum spin type " // char(buffer2))
       end select

    case ("laser")
       input%pertType = iLaser
       call getChildValue(value1, "PolarizationDirection", input%reFieldPolVec)
       call getChildValue(value1, "ImagPolarizationDirection", input%imFieldPolVec, &
            & (/ 0.0_dp, 0.0_dp, 0.0_dp /))
       call getChildValue(value1, "LaserEnergy", input%omega, &
            & modifier=modifier, child=child)
       call convertByMul(char(modifier), energyUnits, child, input%omega)
       call getChildValue(value1, "Phase", input%phase, 0.0_dp)
       call getChildValue(value1, "ExcitedAtoms", buffer, "1:-1", child=child, &
            &multiple=.true.)
       call convAtomRangeToInt(char(buffer), geo%speciesNames, geo%species, &
            &child, input%indExcitedAtom)

       input%nExcitedAtom = size(input%indExcitedAtom)
       if (input%nExcitedAtom == 0) then
          call error("No atoms specified for laser excitation.")
       end if

    case ("kickandlaser")
       input%pertType = iKickAndLaser
       call getChildValue(value1, "KickPolDir", input%polDir)
       if (input%polDir > 4) then
          call detailedError(child, "Wrong specified polarization direction")
       end if
       call getChildValue(value1, "SpinType", input%spType, iTDSinglet)

       call getChildValue(value1, "LaserPolDir", input%reFieldPolVec)
       call getChildValue(value1, "LaserImagPolDir", input%imFieldPolVec, &
            & (/ 0.0_dp, 0.0_dp, 0.0_dp /))
       call getChildValue(value1, "LaserEnergy", input%omega, &
            & modifier=modifier, child=child)
       call convertByMul(char(modifier), energyUnits, child, input%omega)
       call getChildValue(value1, "Phase", input%phase, 0.0_dp)
       call getChildValue(value1, "LaserStrength", input%tdLaserField, modifier=modifier,&
           & child=child)
       call convertByMul(char(modifier), EFieldUnits, child, input%tdLaserField)

       call getChildValue(value1, "ExcitedAtoms", buffer, "1:-1", child=child, &
            &multiple=.true.)
       call convAtomRangeToInt(char(buffer), geo%speciesNames, geo%species, &
            &child, input%indExcitedAtom)
       input%nExcitedAtom = size(input%indExcitedAtom)
       if (input%nExcitedAtom == 0) then
          call error("No atoms specified for laser excitation.")
       end if

    case ("none")
       input%pertType = iNoTDPert

    case default
       call detailedError(child, "Unknown perturbation type " // char(buffer))
    end select

    !! Different envelope functions
    call getChildValue(node, "EnvelopeShape", value1, "Constant")
    call getNodeName(value1, buffer)
    select case(char(buffer))

    case("constant")
       input%envType = iTDConstant

    case("gaussian")
       input%envType = iTDGaussian
       call getChildValue(value1, "Time0", input%time0, 0.0_dp, modifier=modifier, child=child)
       call convertByMul(char(modifier), timeUnits, child, input%Time0)

       call getChildValue(value1, "Time1", input%time1, modifier=modifier, child=child)
       call convertByMul(char(modifier), timeUnits, child, input%Time1)

    case("sin2")
       input%envType = iTDSin2
       call getChildValue(value1, "Time0", input%time0, 0.0_dp, modifier=modifier, child=child)
       call convertByMul(char(modifier), timeUnits, child, input%Time0)

       call getChildValue(value1, "Time1", input%time1, modifier=modifier, child=child)
       call convertByMul(char(modifier), timeUnits, child, input%Time1)

    case("fromfile")
       input%envType = iTDFromFile
       call getChildValue(value1, "Time0", input%time0, 0.0_dp, modifier=modifier, child=child)
       call convertByMul(char(modifier), timeUnits, child, input%Time0)

    case default
       call detailedError(value1, "Unknown envelope shape " // char(buffer))
    end select

    !! Non-adiabatic molecular dynamics
    call getChildValue(node, "IonDynamics", input%tIons, .false.)
    if (input%tIons) then
       call getChildValue(node, "MovedAtoms", buffer, "1:-1", child=child, &
            &multiple=.true.)
       call convAtomRangeToInt(char(buffer), geo%speciesNames, geo%species, &
            &child, input%indMovedAtom)

       input%nMovedAtom = size(input%indMovedAtom)
!       if (input%nMovedAtom == 0) then
!          call error("No atoms specified for molecular dynamics.")
!       end if
       call readInitialVelocitiesNAMD(node, input, geo%nAtom)
       if (input%tReadMDVelocities) then
          ! without a thermostat, if we know the initial
          ! velocities, we do not need a temperature, so just set it to something
          ! 'safe'
          input%tempAtom = minTemp
       else
          call getChildValue(node, "InitialTemperature", input%tempAtom, &
               &modifier=modifier, child=child)
          if (input%tempAtom < 0.0_dp) then
             call detailedError(child, "Negative temperature")
          end if
          call convertByMul(char(modifier), energyUnits, child, input%tempAtom)
          if (input%tempAtom < minTemp) then
             input%tempAtom = 0.0_dp !previously it was minTemp
          end if
       end if
       call getInputMasses(node, geo, masses)
    end if

  end subroutine readElecDynamics


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
    type(listRealR1) :: realBuffer
    integer          :: nVelocities
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
          call convertByMul(char(modifier), VelocityUnits, child, &
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

    type(fnode), pointer :: pGeom, pDevice, pNode, pTask, pTaskType
    type(string) :: buffer, modif
    type(fnode), pointer :: pTmp, field
    type(fnodeList), pointer :: pNodeList
    integer :: ii, contact
    real(dp) :: acc, contactRange(2), lateralContactSeparation
    type(listInt) :: li

    transpar%defined = .true.
    transpar%tPeriodic1D = .not. geom%tPeriodic
    call getChild(root, "Device", pDevice)
    call getChildValue(pDevice, "AtomRange", transpar%idxdevice)
    call getChild(pDevice, "FirstLayerAtoms", pTmp, requested=.false.)
    call readFirstLayerAtoms(pTmp, transpar%PL, transpar%nPLs, transpar%idxdevice)
    if (.not.associated(pTmp)) then
      call setChildValue(pDevice, "FirstLayerAtoms", transpar%PL)
    end if

    call getChild(pDevice, "ContactPLs", pTmp, requested=.false.)
    if (associated(pTmp)) then
      call init(li)
      call getChildValue(pTmp, "", li)
      allocate(transpar%cblk(len(li)))
      call asArray(li,transpar%cblk)
      call destruct(li)
    end if

    !! Note: we parse first the task because we need to know it to defined the
    !! mandatory contact entries. On the other hand we need to wait that
    !! contacts are parsed to resolve the name of the contact for task =
    !! contacthamiltonian
    call getChildValue(root, "Task", pTaskType, child=pTask, default='uploadcontacts')
    call getNodeName(pTaskType, buffer)

    call getChildren(root, "Contact", pNodeList)
    transpar%ncont = getLength(pNodeList)
    if (transpar%ncont < 2) then
      call detailedError(root, "At least two contacts must be defined")
    end if
    allocate(transpar%contacts(transpar%ncont))
    !! Parse contact geometry

    call readContacts(pNodeList, transpar%contacts, geom, (buffer .eq. "uploadcontacts"))

    select case (char(buffer))

    case ("contacthamiltonian")

      transpar%taskUpload = .false.
      call getChildValue(pTaskType, "ContactId", buffer, child=pTmp)
      contact = getContactByName(transpar%contacts(:)%name, tolower(trim(unquote(char(buffer)))),&
          & pTmp)
      transpar%taskContInd = contact
      transpar%contacts(contact)%output = "shiftcont_" // trim(transpar%contacts(contact)%name) //&
          & ".dat"
      if (.not. geom%tPeriodic) then
        call getChildValue(pTaskType, "ContactSeparation", lateralContactSeparation, 1000.0_dp,&
            & modifier=modif, child=field)
        call convertByMul(char(modif),lengthUnits,field,lateralContactSeparation)
      end if
      transpar%tPeriodic1D = .not. geom%tPeriodic

      call reduceGeometry(transpar%contacts(contact)%lattice, transpar%contacts(contact)%idxrange,&
          & lateralContactSeparation, geom)

      transpar%ncont = 0

    case ("uploadcontacts")

      transpar%taskUpload = .true.

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
    integer :: ind, ii

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
      newLatVecs(modulo(ind+1,3)+1, 2) = -newLatVecs(ind,1)
      newLatVecs(ind,2) = 0.0_dp !newLatVecs(modulo(ind+1,3)+1, 1)
      newLatVecs(modulo(ind-1,3)+1, 2) = 0.0_dp
      call cross3(newLatVecs(:,3), newLatVecs(:,1), newLatVecs(:,2))
      newLatVecs(:,2) = newLatVecs(:,2) / sqrt(sum(newLatVecs(:,2)**2))
      newLatVecs(:,3) = newLatVecs(:,3) / sqrt(sum(newLatVecs(:,3)**2))
      newOrigin = 0.0_dp
    end if
    call reduce(geom, contactRange(1), contactRange(2))
    if (.not. geom%tPeriodic) then
      do ii = 2, 3
        minProj = 0_dp !minval(matmul(newLatVecs(:,ii), geom%coords))
        maxProj = 0_dp !maxval(matmul(newLatVecs(:,ii), geom%coords))
        newLatVecs(:,ii) = ((maxProj - minProj) + lateralContactSeparation) * newLatVecs(:,ii)
      end do
    end if
    call setLattice(geom, newOrigin, newLatVecs)

  end subroutine reduceGeometry


  subroutine readFirstLayerAtoms(pnode, pls, npl, idxdevice, check)
    logical, optional :: check
    type(fnode), pointer, intent(in) :: pnode
    integer :: idxdevice(2)
    integer, allocatable :: pls(:)
    integer :: npl

    type(listInt) :: li
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
               &// i2c(idxdevice(1)) // " &
               & and " // i2c(idxdevice(2)) // ".")
          end if
        end if
      else
         npl = 1
         allocate(pls(npl))
         pls = (/ 1 /)
      end if

  end subroutine readFirstLayerAtoms


  subroutine readGreensFunction(pNode, greendens, transpar, tempElec)
    type(TNEGFGreenDensInfo), intent(inout) :: greendens
    type(TTransPar), intent(inout) :: transpar
    real(dp), intent(in) :: tempElec

    type(fnode), pointer :: pGeom, pDevice, pNode, pTask, pTaskType
    type(fnodeList), pointer :: pNodeList
    type(fnode), pointer :: pTmp, field, child1, child2
    real(dp) :: Estep
    integer :: defValue, ii, nfermi
    type(string) :: buffer, modif
    logical :: realAxisConv, equilibrium

    type(listInt) :: li
    type(listReal) :: fermiBuffer

    greendens%defined = .true.

    if (.not. transpar%defined) then
      !! Fermi level: in case of collinear spin we accept two values
      !! (up and down)
      call init(fermiBuffer)
      call getChildValue(pNode, "FermiLevel", fermiBuffer, modifier=modif)
      if ( len(fermiBuffer) .eq. 1) then
        call asArray(fermiBuffer, greendens%oneFermi)
        greendens%oneFermi(2) = greendens%oneFermi(1)
      else if ( len(fermiBuffer) .eq. 2) then
        call asArray(fermiBuffer, greendens%oneFermi)
      else
        call detailedError(pNode, &
            & "FermiLevel accepts 1 or 2 (for collinear spin) values")
      end if
      call destruct(fermiBuffer)
      call convertByMul(char(modif), energyUnits, pNode, greendens%oneFermi)

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
    call getChildValue(pNode, "Delta", greendens%delta, 1.0e-5_dp, modifier=modif, child=field)
    call convertByMul(char(modif), energyUnits, field, greendens%delta)
    call getChildValue(pNode, "SaveSurfaceGFs", greendens%saveSGF, .true.)
    call getChildValue(pNode, "ReadSurfaceGFs", greendens%readSGF, .false.)
    call getChildValue(pNode, "ContourPoints", greendens%nP(1:2), [ 20, 20 ])
    call getChildValue(pNode, "EnclosedPoles",  greendens%nPoles, 3)
    call getChildValue(pNode, "LowestEnergy", greendens%enLow, -2.0_dp, modifier=modif, child=field)
    call convertByMul(char(modif), energyUnits, field, greendens%enLow)
    call getChildValue(pNode, "FermiCutoff", greendens%nkT, 10)
      ! Fermi energy had not been set by other means yet

      ! Non equilibrium integration along real axis:
      ! The code will perform the integration if the number of points is larger
      ! than zero, no matter if there's bias or not.
      ! Therefore I restored the default on the energy step, as it works at zero
      ! bias and it scales flawlessy with increasing bias
      ! It is still allowed to directly set the number of points, if prefered
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
             & modifier=modif)
        call convertByMul(char(modif), energyUnits, child2, Estep)
        realAxisConv = .true.
      ! If the system is under equilibrium we set the number of
      ! points to zero
      else if (equilibrium) then
        call getChildValue(pNode, "RealAxisPoints", greendens%nP(3), &
          & 0, child=child1)
      else
        !Default is a point every 1500H
        call getChildValue(pNode, "RealAxisStep", Estep, 6.65e-4_dp, &
                          &modifier=modif, child=child2)
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


  !> Read in Poisson related data
  subroutine readPoisson(pNode, poisson, tPeriodic, tPeriodic1D)
    type(fnode), pointer :: pNode
    type(TPoissonInfo), intent(inout) :: poisson
    logical, intent(in) :: tPeriodic, tPeriodic1D

    type(fnode), pointer :: pTmp, pTmp2, pChild, field
    type(string) :: buffer, modif
    character(lc) :: strTmp
    real(dp) :: denstol, gatelength_l
    integer :: ii, ibc, bctype
    logical :: needsPoissonBox

    poisson%defined = .true.
    needsPoissonBox = (.not. tPeriodic) .or. tPeriodic1D
    if (needsPoissonBox) then
      call getChildValue(pNode, "PoissonBox", poisson%poissBox, &
          & modifier=modif, child=field)
      call convertByMul(char(modif), lengthUnits, field, &
          & poisson%poissBox)
    end if
    poisson%foundBox = needsPoissonBox
    call getChildValue(pNode, "MinimalGrid", poisson%poissGrid, &
        & [ 0.3_dp, 0.3_dp, 0.3_dp ], modifier=modif, child=field)
    call convertByMul(char(modif), lengthUnits, field, &
        & poisson%poissGrid)
    call getChildValue(pNode, "NumericalNorm", poisson%numericNorm, .false.)
    call getChild(pNode, "AtomDensityCutoff", pTmp, requested=.false., &
        & modifier=modif)
    call getChild(pNode, "AtomDensityTolerance", pTmp2, requested=.false.)
    if (associated(pTmp) .and. associated(pTmp2)) then
      call detailedError(pNode, "Either one of the tags AtomDensityCutoff or&
          & AtomDensityTolerance can be specified.")
    else if (associated(pTmp)) then
      call getChildValue(pTmp, "", poisson%maxRadAtomDens, default=14.0_dp, &
          &  modifier=modif)
      call convertByMul(char(modif), lengthUnits, pTmp, poisson%maxRadAtomDens)
      if (poisson%maxRadAtomDens <= 0.0_dp) then
        call detailedError(pTmp2, "Atom density cutoff must be > 0")
      end if
    else
      call getChildValue(pNode, "AtomDensityTolerance", denstol, 1e-6_dp, &
          & child=pTmp2)
      if (denstol <= 0.0_dp) then
        call detailedError(pTmp2, "Atom density tolerance must be > 0")
      end if
      ! Negative value to signalize automatic determination
      poisson%maxRadAtomDens = -denstol
    end if

    call getChildValue(pNode, "CutoffCheck", poisson%cutoffcheck,&
        & .true.)
    call getChildValue(pNode, "Verbosity", poisson%verbose, 51)
    call getChildValue(pNode, "SavePotential", poisson%savePotential,&
        & .false.)
    call getChildValue(pNode, "PoissonAccuracy", poisson%poissAcc,&
        & 1.0e-6_dp)
    call getChildValue(pNode, "BuildBulkPotential", poisson%bulkBC,&
        & .true.)
    call getChildValue(pNode, "ReadOldBulkPotential", &
        & poisson%readBulkPot, .false.)
    call getChildValue(pNode, "RecomputeAfterDensity",&
        & poisson%solvetwice, .false.)
    call getChildValue(pNode, "MaxPoissonIterations",&
        & poisson%maxPoissIter, 60)

    call getChild(pNode, "OverrideDefaultBC", pTmp, requested=.false.)
    poisson%overrideBC(:) = 0
    if (associated(pTmp)) then
      call getPoissonBoundaryConditionOverrides(pTmp, [ 1, 2 ], &
          & poisson%overrideBC)
    end if

    call getChildValue(pNode, "OverrideBulkBC", pTmp, "none")
    poisson%overrBulkBC(:) = -1
    if (associated(pNode)) then
      call getPoissonBoundaryConditionOverrides(pTmp, [ 0, 1, 2 ], &
          & poisson%overrBulkBC)
    end if

    call getChildValue(pNode, "BoundaryRegion", pTmp, "global")
    call getNodeName(pTmp, buffer)
    select case(char(buffer))
    case ("global")
      poisson%localBCType = "G"
    case ("square")
      poisson%localBCType = "S"
      call getChildValue(pTmp, "BufferLength", poisson%bufferLocBC, &
          &9.0_dp, modifier=modif, child=field)
      call convertByMul(char(modif), lengthUnits, field, &
          & poisson%bufferLocBC)
    case ("circle")
      poisson%localBCType = "C"
      call getChildValue(pTmp, "BufferLength", poisson%bufferLocBC, &
          &9.0_dp, modifier=modif, child=field)
      call convertByMul(char(modif), lengthUnits, field, poisson%bufferLocBC)
    case default
      call getNodeHSDName(pTmp, buffer)
      call detailedError(pTmp, "Invalid boundary region type '" &
          &// char(buffer) // "'")
    end select

    call getChildValue(pNode, "BoxExtension", poisson%bufferBox, &
         &0.0_dp, modifier=modif, child=field)
    call convertByMul(char(modif), lengthUnits, field, poisson%bufferBox)
    if (poisson%bufferBox.lt.0.0_dp) then
      call detailedError(pNode, "BoxExtension must be a positive number")
    endif

    ! PARSE GATE OPTIONS
    call getChildValue(pNode,"Gate",pTmp2,"none",child=pChild)
    call getNodeName(pTmp2, buffer)

    select case(char(buffer))
    case ("none")
      poisson%gateType = "N"
    case ("planar")
      poisson%gateType = "P"
      call getChildValue(pTmp2, "GateLength", poisson%gateLength_l,&
          & 0.0_dp, modifier= modif, child=field)
      call convertByMul(char(modif), lengthUnits, field, &
          &poisson%gateLength_l)

      gatelength_l = poisson%gateLength_l !avoids a warning on intents
      call getChildValue(pTmp2, "GateLength_l", poisson%gateLength_l, &
          & gateLength_l, modifier=modif, child=field)
      call convertByMul(char(modif), lengthUnits, field, &
          &poisson%gateLength_l)

      call getChildValue(pTmp2, "GateLength_t", poisson%gateLength_t, &
          &poisson%gateLength_l, modifier=modif, child=field)
      call convertByMul(char(modif), lengthUnits, field, &
          &poisson%gateLength_t)

      call getChildValue(pTmp2, "GateDistance", poisson%gateRad, &
          &0.0_dp, modifier=modif, child=field)
      call convertByMul(char(modif), lengthUnits, field, &
          &poisson%gateRad)

      call getChildValue(pTmp2, "GatePotential", poisson%gatepot, &
          &0.0_dp, modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, &
          &poisson%gatepot)

      !call getChildValue(pTmp2, "GateDirection", poisson%gatedir, 2)
      poisson%gatedir = 2

    case ("cylindrical")
      poisson%gateType = "C"
      call getChildValue(pTmp2, "GateLength",poisson%gateLength_l,&
          & 0.0_dp, modifier= modif, child=field)
      call convertByMul(char(modif), lengthUnits, field, &
          &poisson%gateLength_l)

      call getChildValue(pTmp2, "GateRadius", poisson%gateRad, &
          &0.0_dp, modifier=modif, child=field)
      call convertByMul(char(modif), lengthUnits, field, &
          &poisson%gateRad)

      call getChildValue(pTmp2, "GatePotential", poisson%gatepot, &
          &0.0_dp, modifier=modif, child=field)
      call convertByMul(char(modif), lengthUnits, field, &
          &poisson%gatepot)

    case default
      call getNodeHSDName(pTmp2, buffer)
      call detailedError(pTmp2, "Invalid gate type '" &
          &// char(buffer) // "'")

    end select

    call getChildValue(pNode, "MaxParallelNodes", poisson%maxNumNodes, 1)

    poisson%scratch = "contacts"

  end subroutine readPoisson


  subroutine getPoissonBoundaryConditionOverrides(pNode, availableConditions, overrideBC)
    type(fnode), pointer, intent(in) :: pNode
    integer, intent(in) :: availableConditions(:)
    integer, intent(inout) :: overrideBC(:)

    integer, parameter :: PERIODIC_BC = 0
    integer, parameter :: DIRICHLET_BC = 1
    integer, parameter :: NEUMANN_BC = 2
    character(10), parameter :: bcstr(0:2) = &
        & [ character(10) :: "Periodic", "Dirichlet", "Neumann" ]
    integer :: bctype, iBC
    integer :: faceBC, oppositeBC
    integer :: ii
    type(listString) :: lStr
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
    type(fnode), pointer :: node
    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb
    !> Atomic geometry, including the contact atoms
    type(TGeometry), intent(in) :: geom
    !> Parameters of the transport calculation
    type(TTransPar), intent(inout) :: tp
    !> Parameters of tunneling and dos calculation
    type(TNEGFTunDos), intent(inout) :: tundos

    type(string) :: model
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
      call detailedError(dephModel,"unkown model")
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

  !!-----------------------------------------------------------------------------
  !! Read Coupling
  !! 2 modes support, constant or specified per each orbital
  !!-----------------------------------------------------------------------------
  subroutine readCoupling(node, elph, geom, orb, tp)
    type(fnode), pointer :: node
    !> container for buttiker-probes parameters
    type(TElPh), intent(inout) :: elph
    !> Geometry type
    type(TGeometry), intent(in) :: geom
    !> Orbitals infos
    type(TOrbitals), intent(in) :: orb
    !> Transport parameter type
    type(TTransPar), intent(in) :: tp

    type(string) :: buffer, method, modif, modif2
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
        & allowEmptyValue=.true., modifier=modif, dummyValue=.true., list=.false.)

    call getNodeName(val, method)

    ! This reads also things like:  "Coupling [eV] = 0.34"
    !if (is_numeric(char(method))) then
    !  call getChildValue(node, "Coupling", rTmp, child=field)
    !  call convertByMul(char(modif), energyUnits, field, rTmp)
    !  elph%coupling = rTmp
    !  return
    !end if

    select case (char(method))
    case ("allorbitals")
      call getChild(child, "AllOrbitals", child2, requested=.false.)
      call getChildValue(child2, "", elph%coupling, child=field)
      call convertByMul(char(modif), energyUnits, field, elph%coupling)

    case ("atomcoupling")
      call getChild(child, "AtomCoupling", child2, requested=.false.)
      allocate(atmCoupling(atm_range(2)-atm_range(1)+1))
      atmCoupling = 0.d0
      call getChildren(child2, "AtomList", children)
      do ii = 1, getLength(children)
        call getItem1(children, ii, child3)
        call getChildValue(child3, "Atoms", buffer, child=child4, &
            &multiple=.true.)
        call convAtomRangeToInt(char(buffer), geom%speciesNames, &
            &geom%species, child4, tmpI1)
        call getChildValue(child3, "Value", rTmp, child=field, modifier=modif2)
        ! If not defined, use common unit modifier defined after Coupling
        if (len(modif2)==0) then
          call convertByMul(char(modif), energyUnits, field, rTmp)
        else
          call convertByMul(char(modif2), energyUnits, field, rTmp)
        end if
        do jj=1, size(tmpI1)
          iAt = tmpI1(jj)
          if (atmCoupling(iAt) /= 0.0_dp) then
            call detailedWarning(child3, "Previous setting of coupling &
                &for atom" // i2c(iAt) // " has been overwritten")
          end if
          atmCoupling(iAt) = rTmp
        enddo
      enddo
      ! Transform atom coupling in orbital coupling
      norbs = 0
      do ii=atm_range(1), atm_range(2)
        elph%coupling(norbs + 1:norbs + orb%nOrbAtom(ii)) = atmCoupling(ii)
        norbs = norbs + orb%nOrbAtom(ii)
      enddo
      deallocate(atmCoupling)

    case ("constant")
      call getChildValue(child, "Constant", rtmp, child=field)
      call convertByMul(char(modif), energyUnits, field, rTmp)
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

    type(fnode), pointer :: pTmp, field
    type(fnode), pointer :: pGeom, pDevice, pNode
    type(fnodeList), pointer :: pNodeList
    integer :: ii, jj, ind, ncont, nKT
    real(dp) :: eRange(2), eRangeDefault(2)
    type(string) :: buffer, modif
    type(WrappedInt1), allocatable :: iAtInRegion(:)
    logical, allocatable :: tShellResInRegion(:)
    character(lc), allocatable :: regionLabelPrefixes(:)
    type(listReal) :: temperature

    tundos%defined = .true.

    ! ncont is needed for contact option allocation
    ncont = transpar%ncont

    call getChildValue(root, "Verbosity", tundos%verbose, 51)
    call getChildValue(root, "WriteLDOS", tundos%writeLDOS, .true.)
    call getChildValue(root, "WriteTunn", tundos%writeTunn, .true.)

    ! Read Temperature. Can override contact definition
    allocate(tundos%kbT(ncont))
    call getChild(root, "ContactTemperature", pTmp, modifier=modif, requested=.false.)
    if (associated(pTmp)) then
      call init(temperature)
      call getChildValue(pTmp, "", temperature)
      if (len(temperature) .ne. ncont) then
        call detailedError(root, "ContactTemperature does not match the number of contacts")
      end if
      call asArray(temperature, tundos%kbT)
      call destruct(temperature)
      call convertByMul(char(modif), energyUnits, pTmp, tundos%kbT)
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
      call getChildValue(root, "EnergyRange", eRange, modifier=modif,&
      & child=field)
      call convertByMul(char(modif), energyUnits, field, eRange)
      call getChildValue(root, "EnergyStep", tundos%estep,&
      & modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, tundos%estep)
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
                          &modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, tundos%estep)
      call getChildValue(root, "EnergyRange", eRange, eRangeDefault, &
                          modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, eRange)
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
          call getEmitterCollectorByName(pNode, tundos%ni(ii),&
              & tundos%nf(ii), transpar%contacts(:)%name)
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
          &1.0e-5_dp, modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, &
          &tundos%delta)
      call getChildValue(root, "BroadeningDelta", tundos%broadeningDelta, &
          &0.0_dp, modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, &
          &tundos%broadeningDelta)

      call getChildren(root, "Region", pNodeList)
      call readPDOSRegions(pNodeList, geo, iAtInRegion, tShellResInRegion, &
          & regionLabelPrefixes)
      call destroyNodeList(pNodeList)
      call transformPdosRegionInfo(iAtInRegion, tShellResInRegion, &
          & regionLabelPrefixes, orb, geo%species, tundos%dosOrbitals, &
          & tundos%dosLabels)

  end subroutine readTunAndDos

  !> Read bias information, used in Analysis and Green's function eigensolver
  subroutine readContacts(pNodeList, contacts, geom, upload)
    type(ContactInfo), allocatable, dimension(:), intent(inout) :: contacts
    type(fnodeList), pointer :: pNodeList
    type(TGeometry), intent(in) :: geom
    logical, intent(in) :: upload

    real(dp) :: contactLayerTol
    integer :: ii, jj
    type(fnode), pointer :: field, pNode, pTmp, pWide
    type(string) :: buffer, modif
    type(listReal) :: fermiBuffer


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
        call detailedError(pTmp, "Contact id '" // trim(contacts(ii)%name) &
            &//  "' already in use")
      end if

      call getChildValue(pNode, "PLShiftTolerance", contactLayerTol, 1e-5_dp, modifier=modif,&
          & child=field)
      call convertByMul(char(modif), lengthUnits, field, contactLayerTol)
      call getChildValue(pNode, "AtomRange", contacts(ii)%idxrange, child=pTmp)
      call getContactVector(contacts(ii)%idxrange, geom, ii, contacts(ii)%name, pTmp,&
          & contactLayerTol, contacts(ii)%lattice, contacts(ii)%dir)
      contacts(ii)%length = sqrt(sum(contacts(ii)%lattice**2))

      ! Contact temperatures. A negative default is used so it is quite clear when the user sets a
      ! different value. In such a case this overrides values defined in the Filling block
      call getChild(pNode,"Temperature", field, modifier=modif, requested=.false.)
      if (associated(field)) then
        call getChildValue(pNode, "Temperature", contacts(ii)%kbT, 0.0_dp, modifier=modif,&
            & child=field)
        call convertByMul(char(modif), energyUnits, field, contacts(ii)%kbT)
      else
        contacts(ii)%kbT = -1.0_dp ! -1.0 simply means 'not defined'
      end if

      if (upload) then
        call getChildValue(pNode, "Potential", contacts(ii)%potential, 0.0_dp, modifier=modif,&
            & child=field)
        call convertByMul(char(modif), energyUnits, field, contacts(ii)%potential)

        call getChildValue(pNode, "WideBand", contacts(ii)%wideBand, .false.)

        if (contacts(ii)%wideBand) then

          ! WideBandApproximation is defined as energy spacing between levels of the contact. In the
          ! code the inverse value (Density of states) is used. Convert the negf input
          ! value. Default is 20 / e eV.
          call getChildValue(pNode, "LevelSpacing", contacts(ii)%wideBandDos, 0.735_dp,&
              & modifier=modif, child=field)
          call convertByMul(char(modif), energyUnits, field, contacts(ii)%wideBandDos)
          contacts(ii)%wideBandDos = 1.d0 / contacts(ii)%wideBandDos

        end if

        ! Fermi level: in case of collinear spin we accept two values (up and down)
        call init(fermiBuffer)
        call getChildValue(pNode, "FermiLevel", fermiBuffer, modifier=modif)
        if ( len(fermiBuffer) .eq. 1) then
          call asArray(fermiBuffer, contacts(ii)%eFermi)
          contacts(ii)%eFermi(2) = contacts(ii)%eFermi(1)
        else if ( len(fermiBuffer) .eq. 2) then
          call asArray(fermiBuffer, contacts(ii)%eFermi)
        else
          call detailedError(pNode, "FermiLevel accepts 1 or 2 (for collinear spin) values")
        end if
        call destruct(fermiBuffer)
        call convertByMul(char(modif), energyUnits, pNode, contacts(ii)%eFermi)
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

    end do

  end subroutine readContacts



  !> Read in Fermi levels
  subroutine getFermiLevels(pNode, eFermis, nodeModifier)
    type(fnode), pointer :: pNode
    real(dp), intent(out) :: eFermis(:)
    type(string), intent(in) :: nodeModifier

    real(dp) :: eFermi
    type(fnode), pointer :: pChild
    type(string) :: modifier

    call getChild(pNode, "SetForAll", pChild, requested=.false.)
    if (associated(pChild)) then
      call getChildValue(pChild, "", eFermi)
      call convertByMul(char(nodeModifier), energyUnits, pNode, eFermi)
      eFermis(:) = eFermi
    else
      call getChildValue(pNode, "", eFermis, modifier=modifier,&
          & child=pChild)
      call convertByMul(char(modifier), energyUnits, pChild, eFermis)
    end if

  end subroutine getFermiLevels


  !> Get contacts for terminal currents by name
  subroutine getEmitterCollectorByName(pNode, emitter, collector, contactNames)
    type(fnode), pointer :: pNode
    integer, intent(out) :: emitter, collector
    character(len=*), intent(in) :: contactNames(:)

    type(listString) :: lString
    character(len=mc) :: buffer
    integer :: ind
    logical :: tFound

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
    character(len=*), intent(in) :: contactNames(:)
    character(len=*), intent(in) :: contName
    type(fnode), pointer :: pNode
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
  subroutine readPDOSRegions(children, geo, iAtInregion, tShellResInRegion, regionLabels)
    type(fnodeList), pointer :: children
    type(TGeometry), intent(in) :: geo
    type(WrappedInt1), allocatable, intent(out) :: iAtInRegion(:)
    logical, allocatable, intent(out) :: tShellResInRegion(:)
    character(lc), allocatable, intent(out) :: regionLabels(:)

    integer :: nReg, iReg
    integer, allocatable :: tmpI1(:)
    type(fnode), pointer :: child, child2
    type(string) :: buffer
    character(lc) :: strTmp

    nReg = getLength(children)
    allocate(tShellResInRegion(nReg))
    allocate(regionLabels(nReg))
    allocate(iAtInRegion(nReg))
    do iReg = 1, nReg
      call getItem1(children, iReg, child)
      call getChildValue(child, "Atoms", buffer, child=child2, &
          & multiple=.true.)
      call convAtomRangeToInt(char(buffer), geo%speciesNames, &
          & geo%species, child2, tmpI1)
      iAtInRegion(iReg)%data = tmpI1
      call getChildValue(child, "ShellResolved", &
          & tShellResInRegion(iReg), .false., child=child2)
      if (tShellResInRegion(iReg)) then
        if (.not. all(geo%species(tmpI1) == geo%species(tmpI1(1)))) then
          call detailedError(child2, "Shell resolved PDOS can only summed up &
              & over atoms of the same type")
        end if
      end if
      write(strTmp, "('region',I0)") iReg
      call getChildValue(child, "Label", buffer, trim(strTmp))
      regionLabels(iReg) = unquote(char(buffer))
    end do

  end subroutine readPDOSRegions


  !> Some assignment and consistency check in negf/poisson containers before calling initialization
  subroutine finalizeNegf(input)
    type(inputData), intent(inout) :: input

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

    !! Temporarily not supporting surface green function read/load
    !! for spin polarized, because spin is handled outside of libnegf
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


  !> Reads the parallel block.
  subroutine readParallel(root, input)

    !> Root node eventually containing the current block
    type(fnode), pointer, intent(in) :: root

    !> Input structure to be filled
    type(inputData), intent(inout) :: input

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
    #:if WITH_TRANSPORT
      if (input%transpar%ncont > 0 .and. input%ctrl%parallelOpts%nGroup > 1) then
        call detailedError(pTmp, "Multiple processor groups are currently incompatible with&
            & transport.")
      end if
    #:endif
      call getChildValue(node, "UseOmpThreads", input%ctrl%parallelOpts%tOmpThreads, .not. withMpi)
      call readBlacs(node, input%ctrl%parallelOpts%blacsOpts)
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

  function is_numeric(string) result(is)
    character(len=*), intent(in) :: string
    logical :: is

    real :: x
    integer :: err

    print*,string

    read(string,*,iostat=err) x
    is = (err == 0)
    print*, x, err, is

  end function is_numeric

end module parser
