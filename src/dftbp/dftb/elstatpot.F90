!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"
#:include 'error.fypp'

!> Provides data structure for evaluated electrostatic potentials
module dftbp_dftb_elstatpot
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_environment, only : TEnvironment
  use dftbp_common_hamiltoniantypes, only : hamiltonianTypes
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_extfields, only : TEField
  use dftbp_dftb_scc, only : TScc
  implicit none

  private
  public :: TElStatPotentialsInp
  public :: TElStatPotentials, TElStatPotentials_init


  !> Electrostatic potential evaluator
  type :: TElStatPotentialsInp

    !> File to store the resulting points
    character(lc) :: espOutFile = 'ESP.dat'

    !> Should the potential appended to the file
    logical :: tAppendESP = .false.

    !> Location of electrostatic potential points
    real(dp), allocatable :: espGrid(:,:)

    !> Size of the grid if regular, 0 otherwise
    integer :: gridDimensioning(3) = 0

    !> Origin of the grid if regular
    real(dp) :: origin(3)

    !> Axes of the grid if regular
    real(dp) :: axes(3,3)

    !> short range softening of the potential
    real(dp) :: softenEsp = 1.0E-6_dp

  end type TElStatPotentialsInp


  !> Contains the potential
  type :: TElStatPotentials

    !> Points to evaluate the field if requested
    real(dp), allocatable :: espGrid(:,:)

    !> Size of the grid if regular, 0 otherwise
    integer :: gridDimensioning(3) = 0

    !> Origin of the grid if regular
    real(dp) :: origin(3)

    !> Axes of the grid if regular
    real(dp) :: axes(3,3)

    !> Value of a short-distance softening term
    real(dp) :: softenEsp

    !> File containing output potentials
    character(lc) :: espOutFile

    !> should the file be appended or overwritten
    logical :: tAppendEsp

    !> Internal electrostatic potential from DFTB model
    real(dp), allocatable :: intPotential(:)

    !> Contribution from external potential(s)
    real(dp), allocatable :: extPotential(:)

  contains

    !> Calculate potentials
    procedure :: evaluate

  end type TElStatPotentials


contains


  !> Initialises calculator instance.
  subroutine TElStatPotentials_init(this, input, isAnExtPotential, iHamiltonianType, errStatus)

    !> Instance of this calculator
    type(TElStatPotentials), intent(out) :: this

    !> Input data
    type(TElStatPotentialsInp), intent(inout) :: input

    !> Is an external potential being evaluated
    logical, intent(in) :: isAnExtPotential

    !> Choice of electronic hamiltonian
    integer, intent(in) :: iHamiltonianType

    !> Status of operation
    type(TStatus), intent(out) :: errStatus

    select case(iHamiltonianType)
    case default
      @:RAISE_ERROR(errStatus, -1, "Unknown hamiltonian type")
    case(hamiltonianTypes%dftb)

    case(hamiltonianTypes%xtb)
      @:RAISE_ERROR(errStatus, -1, "Electrostatic potential plotting not currently available for&
          & xTB")
    end select

    this%espOutFile = input%espOutFile
    this%tAppendEsp = input%tAppendEsp
    call move_alloc(input%espGrid, this%espGrid)
    this%gridDimensioning = input%gridDimensioning
    this%origin = input%origin
    this%axes = input%axes
    this%softenEsp = input%softenEsp
    allocate(this%intPotential(size(this%espGrid,dim=2)))
    if (isAnExtPotential) then
      allocate(this%extPotential(size(this%espGrid,dim=2)))
    end if


  end subroutine TElStatPotentials_init


  !> Evaluate the electrostatic potential at specified points.
  !>
  !> Note, internally to DFTB+ the potential sign is defined opposite to the usual convention, since
  !> the charge on the electron is not included when applying potentials to the hamiltonian.
  subroutine evaluate(this, env, SccCalc, EField)

    !> Object holding the potential location information
    class(TElStatPotentials), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Module variables for SCC
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Electric field magnitude
    type(TEField), intent(in), optional :: eField

    integer :: ii

    call sccCalc%getInternalElStatPotential(this%intPotential, env, this%espGrid,&
        & epsSoften=this%softenEsp)
    if (allocated(this%extPotential)) then
      call sccCalc%getExternalElStatPotential(this%extPotential, env, this%espGrid,&
          & epsSoften=this%softenEsp)
      this%extPotential = -this%extPotential
      if (present(eField)) then
        if (allocated(eField%EFieldStrength)) then
          do ii = 1, size(this%espGrid,dim=2)
            this%extPotential(ii) = this%extPotential(ii) + dot_product(this%espGrid(:, ii),&
                & eField%eField)
          end do
        end if
      end if
    end if

  end subroutine evaluate


end module dftbp_dftb_elstatpot
