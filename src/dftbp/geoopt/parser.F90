!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Module to read input from HSD tree
module dftbp_geoopt_parser
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_globalenv, only : stdOut
  use dftbp_common_unitconversion, only : timeUnits, lengthUnits, energyUnits, forceUnits
  use dftbp_extlibs_xmlf90, only : fnode, string, char, getNodeName
  use dftbp_geoopt_input, only : TGeoOptInput, TFilterInput, TRationalFunctionInput, &
      & TLBFGSInput, TFireInput, TOptConv
  use dftbp_io_charmanip, only : unquote
  use dftbp_io_hsdutils, only : getChild, getChildValue, setChild, detailedError, &
      & detailedWarning, getSelectedAtomIndices
  use dftbp_io_hsdutils2, only : convertByMul
  use dftbp_type_typegeometry, only : TGeometry
  implicit none
  private
  public :: readGeoOptimizer


contains


  !> General entry point to read geometry optimization
  subroutine readGeoOptimizer(node, input, geom)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TGeoOptInput), intent(out) :: input

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    type(fnode), pointer :: child, value1
    type(string) :: buffer

    call getChildValue(node, "Optimizer", child)
    call getNodeName(child, buffer)

    select case (char(buffer))
    case default
      call detailedError(node, "Invalid optimizer name.")
    case("fire")
      allocate(input%fire)
      call readFire(child, input%fire, geom)
    case("lbfgs")
      allocate(input%lbfgs)
      call readLBFGS(child, input%lbfgs, geom)
    case("rf")
      allocate(input%rf)
      call readRF(child, input%rf, geom)
    end select

    call getChildValue(node, "Filter", child, "Cartesian")
    call getNodeName(child, buffer)

    select case (char(buffer))
    case default
      call detailedError(node, "Invalid filter name.")
    case("cartesian")
      call readCartesianFilter(child, input%filter, geom)
    end select

    call getChildValue(node, "Convergence", value1, "", child=child, allowEmptyValue=.true.)
    call readConvergence(child, input%conv)

    call getChildValue(node, "MaxSteps", input%nGeoSteps, 20*geom%nAtom)

    call getChildValue(node, "OutputPrefix", buffer, "geo_end")
    input%outFile = trim(unquote(char(buffer)))

  end subroutine readGeoOptimizer


  !> Entry point for reading convergence thresholds
  subroutine readConvergence(node, input)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TOptConv), intent(out) :: input

    type(fnode), pointer :: field
    type(string) :: modifier

    call getChildValue(node, "Energy", input%ethr, huge(1.0_dp), &
        & modifier=modifier, child=field)
    call convertByMul(char(modifier), energyUnits, field, input%ethr)

    call getChildValue(node, "GradNorm", input%gthr, huge(1.0_dp), &
        & modifier=modifier, child=field)
    call convertByMul(char(modifier), forceUnits, field, input%gthr)
    call getChildValue(node, "GradAMax", input%gmax, 1.0e-4_dp, &
        & modifier=modifier, child=field)
    call convertByMul(char(modifier), forceUnits, field, input%gmax)

    call getChildValue(node, "DispNorm", input%dthr, huge(1.0_dp), &
        & modifier=modifier, child=field)
    call convertByMul(char(modifier), lengthUnits, field, input%dthr)
    call getChildValue(node, "DispAMax", input%dmax, huge(1.0_dp), &
        & modifier=modifier, child=field)
    call convertByMul(char(modifier), lengthUnits, field, input%dmax)
  end subroutine readConvergence


  !> Entry point for reading input for FIRE
  subroutine readFire(node, input, geom)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TFireInput), intent(out) :: input

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    type(fnode), pointer :: field
    type(string) :: modifier

    call getChildValue(node, "nMin", input%nMin, 5)
    call getChildValue(node, "aPar", input%a_start, 0.1_dp)
    call getChildValue(node, "fInc", input%f_inc, 1.1_dp)
    call getChildValue(node, "fDec", input%f_dec, 0.5_dp)
    call getChildValue(node, "fAlp", input%f_alpha, 0.99_dp)
    call getChildValue(node, "StepSize", input%dt_max, 1.0_dp, &
        &modifier=modifier, child=field)
    call convertByMul(char(modifier), timeUnits, field, input%dt_max)
  end subroutine readFire


  !> Entry point for reading input for LBFGS
  subroutine readLBFGS(node, input, geom)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TLBFGSInput), intent(out) :: input

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    call getChildValue(node, "Memory", input%memory, 20)
  end subroutine readLBFGS


  !> Entry point for reading input for rational function optimizer
  subroutine readRF(node, input, geom)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TRationalFunctionInput), intent(out) :: input

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    call getChildValue(node, "diagLimit", input%diagLimit, 1.0e-2_dp)
  end subroutine readRF


  !> Entry point for reading input for cartesian geometry transformation filter
  subroutine readCartesianFilter(node, input, geom)

    !> Node to get the information from
    type(fnode), pointer :: node

    !> Control structure to be filled
    type(TFilterInput), intent(out) :: input

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    type(fnode), pointer :: child
    type(string) :: buffer

    call getChildValue(node, "LatticeOpt", input%lattice, .false.)
    call getChildValue(node, "MovedAtoms", buffer, "1:-1", multiple=.true., child=child)
    call getSelectedAtomIndices(child, char(buffer), geom%speciesNames, geom%species, &
        & input%indMovedAtom)
  end subroutine readCartesianFilter


end module dftbp_geoopt_parser
