!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'error.fypp'

!> Module to read input from HSD tree
module dftbp_dftbplus_input_geoopt
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_status, only : TStatus
  use dftbp_common_unitconversion, only : timeUnits, lengthUnits, energyUnits, forceUnits
  use dftbp_extlibs_xmlf90, only : fnode, string, char, getNodeName
  use dftbp_geoopt_package, only : TFilterInput, TOptimizerInput, TRationalFuncInput,&
      & TLbfgsInput, TFireInput, TSteepdescInput, TOptTolerance
  use dftbp_io_charmanip, only : unquote
  use dftbp_io_hsdutils, only : getChild, getChildValue, setChild, detailedError, detailedWarning,&
      & getSelectedAtomIndices
  use dftbp_io_hsdutils2, only : convertUnitHsd, renameChildren
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  public :: readGeoOptInput, readOptimizerInput, TGeoOptInput


  !> General input wrapper for optimisers in this package
  type :: TGeoOptInput

    !> Input for coordinate transformation and filter step
    type(TFilterInput) :: filter

    !> Optimiser input choice
    class(TOptimizerInput), allocatable :: optimiser

    !> Tolerances for optimization
    type(TOptTolerance) :: tolerance

    !> Number of allowed geometry optimization steps
    integer :: nGeoSteps = huge(1) - 1

    !> Prefix of the output file name
    character(len=:), allocatable :: outFile

  end type TGeoOptInput


contains

  !> General entry point to read geometry optimization
  subroutine readGeoOptInput(node, geom, input, atomsRange, errStatus)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Control structure to be filled
    type(TGeoOptInput), intent(out) :: input

    !> Default range of moving atoms (may be restricted for example by contacts in transport
    !> calculations)
    character(len=*), intent(in) :: atomsRange

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child, value1
    type(string) :: buffer

    call renameChildren(node, "Optimizer", "Optimiser")
    call getChildValue(node, "Optimiser", child, errStatus, "Rational")
    @:PROPAGATE_ERROR(errStatus)
    call readOptimizerInput(child, input%optimiser, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call readFilterInput(node, geom, input%filter, atomsRange, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(node, "Convergence", value1, errStatus, "", child=child,&
        & allowEmptyValue=.true.)
    @:PROPAGATE_ERROR(errStatus)
    call readOptTolerance(child, input%tolerance, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(node, "MaxSteps", input%nGeoSteps, errStatus, 20*geom%nAtom)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(node, "OutputPrefix", buffer, errStatus, "geo_end")
    @:PROPAGATE_ERROR(errStatus)
    input%outFile = trim(unquote(char(buffer)))

  end subroutine readGeoOptInput


  !> Reads the optimiser
  subroutine readOptimizerInput(node, input, errStatus)

    !> Optimiser node
    type(fnode), pointer, intent(in) :: node

    !> Control structure to be filled
    class(TOptimizerInput), allocatable, intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(TSteepDescInput), allocatable :: steepDescInput
    type(TFireInput), allocatable :: fireInput
    type(TLbfgsInput), allocatable :: lbfgsInput
    type(TRationalFuncInput), allocatable :: rationalFuncInput
    type(string) :: buffer

    call getNodeName(node, buffer)
    select case (char(buffer))
    case default
      call detailedError(node, "Invalid optimiser name.", errStatus)
      @:PROPAGATE_ERROR(errStatus)
    case("steepestdescent")
        allocate(steepDescInput)
        call readSteepDescInput(node, steepDescInput, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call move_alloc(steepDescInput, input)
    case("fire")
        allocate(fireInput)
        call readFireInput(node, fireInput, errStatus)
        @:PROPAGATE_ERROR(errStatus)
        call move_alloc(fireInput, input)
    case("lbfgs")
      allocate(lbfgsInput)
      call readLbfgsInput(node, lbfgsInput, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call move_alloc(lbfgsInput, input)
    case("rational")
      allocate(rationalFuncInput)
      call readRationalFuncInput(node, rationalFuncInput, errStatus)
      @:PROPAGATE_ERROR(errStatus)
      call move_alloc(rationalFuncInput, input)
    end select

  end subroutine readOptimizerInput


  !> Entry point for reading input for cartesian geometry transformation filter
  subroutine readFilterInput(node, geom, input, atomsRange, errStatus)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Control structure to be filled
    type(TFilterInput), intent(out) :: input

    !> Default range of moving atoms (may be restricted for example by contacts in transport
    !> calculations)
    character(len=*), intent(in) :: atomsRange

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: child
    type(string) :: buffer

    call getChildValue(node, "LatticeOpt", input%lattice, errStatus, .false.)
    @:PROPAGATE_ERROR(errStatus)
    if (input%lattice) then
      call getChildValue(node, "FixAngles", input%fixAngles, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
      call getChildValue(node, "FixLengths", input%fixLength, errStatus,&
          & [.false., .false., .false.])
      @:PROPAGATE_ERROR(errStatus)
      if (input%fixAngles .and. all(input%fixLength)) then
        call detailedError(node, "LatticeOpt with all lattice vectors fixed is not possible",&
            & errStatus)
        @:PROPAGATE_ERROR(errStatus)
      end if
      call getChildValue(node, "Isotropic", input%isotropic, errStatus, .false.)
      @:PROPAGATE_ERROR(errStatus)
    end if
    call getChildValue(node, "MovedAtoms", buffer, errStatus, trim(atomsRange), multiple=.true.,&
        & child=child)
    @:PROPAGATE_ERROR(errStatus)
    call getSelectedAtomIndices(child, char(buffer), geom%speciesNames, geom%species,&
        & input%indMovedAtom, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readFilterInput


  !> Entry point for reading convergence thresholds
  subroutine readOptTolerance(node, input, errStatus)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TOptTolerance), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: field
    type(string) :: modifier

    call getChildValue(node, "Energy", input%energy, errStatus, huge(1.0_dp), modifier=modifier,&
        & child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), energyUnits, field, input%energy, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(node, "GradNorm", input%gradNorm, errStatus, huge(1.0_dp),&
        & modifier=modifier, child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), forceUnits, field, input%gradNorm, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "GradElem", input%gradElem, errStatus, 1.0e-4_dp, modifier=modifier,&
        & child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), forceUnits, field, input%gradElem, errStatus)
    @:PROPAGATE_ERROR(errStatus)

    call getChildValue(node, "DispNorm", input%dispNorm, errStatus, huge(1.0_dp),&
        & modifier=modifier, child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), lengthUnits, field, input%dispNorm, errStatus)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "DispElem", input%dispElem, errStatus, huge(1.0_dp),&
        & modifier=modifier, child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), lengthUnits, field, input%dispElem, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readOptTolerance


  !> Entry point for reading input for SteepestDescent
  subroutine readSteepDescInput(node, input, errStatus)

    !> Node to get the information from
    type(fnode), intent(in), pointer :: node

    !> Control structure to be filled
    type(TSteepDescInput), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    call getChildValue(node, "ScalingFactor", input%scalingFactor, errStatus, 1.0_dp)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readSteepDescInput


  !> Entry point for reading input for FIRE
  subroutine readFireInput(node, input, errStatus)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TFireInput), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(fnode), pointer :: field
    type(string) :: modifier

    call getChildValue(node, "nMin", input%nMin, errStatus, 5)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "aPar", input%a_start, errStatus, 0.1_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "fInc", input%f_inc, errStatus, 1.1_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "fDec", input%f_dec, errStatus, 0.5_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "fAlpha", input%f_alpha, errStatus, 0.99_dp)
    @:PROPAGATE_ERROR(errStatus)
    call getChildValue(node, "StepSize", input%dt_max, errStatus, 1.0_dp, modifier=modifier,&
        & child=field)
    @:PROPAGATE_ERROR(errStatus)
    call convertUnitHsd(char(modifier), timeUnits, field, input%dt_max, errStatus)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readFireInput


  !> Entry point for reading input for LBFGS
  subroutine readLbfgsInput(node, input, errStatus)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TLbfgsInput), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    call getChildValue(node, "Memory", input%memory, errStatus, 20)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readLbfgsInput


  !> Entry point for reading input for rational function optimiser
  subroutine readRationalFuncInput(node, input, errStatus)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TRationalFuncInput), intent(out) :: input

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    call getChildValue(node, "DiagLimit", input%diagLimit, errStatus, 1.0e-2_dp)
    @:PROPAGATE_ERROR(errStatus)

  end subroutine readRationalFuncInput


end module dftbp_dftbplus_input_geoopt
