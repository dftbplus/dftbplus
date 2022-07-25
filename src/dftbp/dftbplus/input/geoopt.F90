!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Module to read input from HSD tree
module dftbp_dftbplus_input_geoopt
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_unitconversion, only : timeUnits, lengthUnits, energyUnits, forceUnits
  use dftbp_extlibs_xmlf90, only : fnode, string, char, getNodeName
  use dftbp_geoopt_package, only : TFilterInput, TOptimizerInput, TRationalFuncInput,&
      & TLbfgsInput, TFireInput, TOptTolerance
  use dftbp_io_charmanip, only : unquote
  use dftbp_io_hsdutils, only : getChild, getChildValue, setChild, detailedError, &
      & detailedWarning, getSelectedAtomIndices
  use dftbp_io_hsdutils2, only : convertUnitHsd
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  public :: readGeoOptInput, TGeoOptInput


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
  subroutine readGeoOptInput(node, geom, input, atomsRange)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Control structure to be filled
    type(TGeoOptInput), intent(out) :: input

    !> Default range of moving atoms (may be restricted for example by contacts in transport
    !> calculations)
    character(len=*), intent(in) :: atomsRange

    type(fnode), pointer :: child, value1
    type(string) :: buffer

    call getChildValue(node, "Optimiser", child, "Rational")
    call readOptimizerInput(child, input%optimiser)

    call readFilterInput(node, geom, input%filter, atomsRange)

    call getChildValue(node, "Convergence", value1, "", child=child, allowEmptyValue=.true.)
    call readOptTolerance(child, input%tolerance)

    call getChildValue(node, "MaxSteps", input%nGeoSteps, 20*geom%nAtom)

    call getChildValue(node, "OutputPrefix", buffer, "geo_end")
    input%outFile = trim(unquote(char(buffer)))

  end subroutine readGeoOptInput


  !> Reads the optimiser
  subroutine readOptimizerInput(node, input)

    !> Optimiser node
    type(fnode), pointer, intent(in) :: node

    !> Control structure to be filled
    class(TOptimizerInput), allocatable, intent(out) :: input

    type(TFireInput), allocatable :: fireInput
    type(TLbfgsInput), allocatable :: lbfgsInput
    type(TRationalFuncInput), allocatable :: rationalFuncInput
    type(string) :: buffer

    call getNodeName(node, buffer)
    select case (char(buffer))
    case default
      call detailedError(node, "Invalid optimiser name.")
    case("fire")
        allocate(fireInput)
        call readFireInput(node, fireInput)
        call move_alloc(fireInput, input)
    case("lbfgs")
      allocate(lbfgsInput)
      call readLbfgsInput(node, lbfgsInput)
      call move_alloc(lbfgsInput, input)
    case("rational")
      allocate(rationalFuncInput)
      call readRationalFuncInput(node, rationalFuncInput)
      call move_alloc(rationalFuncInput, input)
    end select

  end subroutine readOptimizerInput


  !> Entry point for reading input for cartesian geometry transformation filter
  subroutine readFilterInput(node, geom, input, atomsRange)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    !> Control structure to be filled
    type(TFilterInput), intent(out) :: input

    !> Default range of moving atoms (may be restricted for example by contacts in transport
    !> calculations)
    character(len=*), intent(in) :: atomsRange

    type(fnode), pointer :: child
    type(string) :: buffer

    call getChildValue(node, "LatticeOpt", input%lattice, .false.)
    if (input%lattice) then
      call getChildValue(node, "FixAngles", input%fixAngles, .false.)
      call getChildValue(node, "FixLengths", input%fixLength, [.false., .false., .false.])
      if (input%fixAngles .and. all(input%fixLength)) then
        call detailedError(node, "LatticeOpt with all lattice vectors fixed is not possible")
      end if
      call getChildValue(node, "Isotropic", input%isotropic, .false.)
    end if
    call getChildValue(node, "MovedAtoms", buffer, trim(atomsRange), multiple=.true., child=child)
    call getSelectedAtomIndices(child, char(buffer), geom%speciesNames, geom%species, &
        & input%indMovedAtom)

  end subroutine readFilterInput


  !> Entry point for reading convergence thresholds
  subroutine readOptTolerance(node, input)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TOptTolerance), intent(out) :: input

    type(fnode), pointer :: field
    type(string) :: modifier

    call getChildValue(node, "Energy", input%energy, huge(1.0_dp), modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), energyUnits, field, input%energy)

    call getChildValue(node, "GradNorm", input%gradNorm, huge(1.0_dp), modifier=modifier,&
        & child=field)
    call convertUnitHsd(char(modifier), forceUnits, field, input%gradNorm)
    call getChildValue(node, "GradElem", input%gradElem, 1.0e-4_dp, modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), forceUnits, field, input%gradElem)

    call getChildValue(node, "DispNorm", input%dispNorm, huge(1.0_dp), modifier=modifier,&
        & child=field)
    call convertUnitHsd(char(modifier), lengthUnits, field, input%dispNorm)
    call getChildValue(node, "DispElem", input%dispElem, huge(1.0_dp), modifier=modifier,&
        & child=field)
    call convertUnitHsd(char(modifier), lengthUnits, field, input%dispElem)

  end subroutine readOptTolerance


  !> Entry point for reading input for FIRE
  subroutine readFireInput(node, input)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TFireInput), intent(out) :: input

    type(fnode), pointer :: field
    type(string) :: modifier

    call getChildValue(node, "nMin", input%nMin, 5)
    call getChildValue(node, "aPar", input%a_start, 0.1_dp)
    call getChildValue(node, "fInc", input%f_inc, 1.1_dp)
    call getChildValue(node, "fDec", input%f_dec, 0.5_dp)
    call getChildValue(node, "fAlpha", input%f_alpha, 0.99_dp)
    call getChildValue(node, "StepSize", input%dt_max, 1.0_dp, modifier=modifier, child=field)
    call convertUnitHsd(char(modifier), timeUnits, field, input%dt_max)

  end subroutine readFireInput


  !> Entry point for reading input for LBFGS
  subroutine readLbfgsInput(node, input)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TLbfgsInput), intent(out) :: input

    call getChildValue(node, "Memory", input%memory, 20)

  end subroutine readLbfgsInput


  !> Entry point for reading input for rational function optimiser
  subroutine readRationalFuncInput(node, input)

    !> Node to get the information from
    type(fnode), pointer, intent(in) :: node

    !> Control structure to be filled
    type(TRationalFuncInput), intent(out) :: input

    call getChildValue(node, "DiagLimit", input%diagLimit, 1.0e-2_dp)

  end subroutine readRationalFuncInput


end module dftbp_dftbplus_input_geoopt
