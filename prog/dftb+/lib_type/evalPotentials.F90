!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "common.fypp"

!> Provides data structure for evaluated electrostatic potentials
module electrostaticPotentials
  use accuracy
  use scc
  use environment
  implicit none
  private

  public :: TElectrostaticPotentialsInp, TElectrostaticPotentials, initialise

  type :: TElectrostaticPotentialsInp

    !> File to store the resulting points
    character(lc) :: EspOutFile = 'ESP.dat'

    !> Should the potential appended to the file
    logical :: tAppendESP = .false.

    !> Location of electrostatic potential points
    real(dp), allocatable :: ESPgrid(:,:)

    !> Size of the grid if regular, 0 otherwise
    integer :: gridDimensioning(3) = 0

    !> short range softening of the potential
    real(dp) :: softenESP = 1.0E-6_dp

  end type TElectrostaticPotentialsInp

  !> Contains the potential
  type :: TElectrostaticPotentials

    !> Points to evaluate the field if requested
    real(dp), allocatable :: ESPgrid(:,:)

    !> Size of the grid if regular, 0 otherwise
    integer :: gridDimensioning(3) = 0

    !> Value of a short-distance softening term
    real(dp) :: softenESP

    !> file unit for ESP result
    integer :: fdESP

    !> File containing output potentials
    character(lc) :: EspOutFile

    !> should the file be appended or overwritten
    logical :: tAppendESP

    real(dp), allocatable :: ESPpotential(:)

    real(dp), allocatable :: extESPpotential(:)


  contains

    !> Calculate potentials
    procedure :: evaluate

  end type TElectrostaticPotentials

  interface initialise
    module procedure TElectrostaticPotentials_initialise
  end interface initialise

contains

  !> Initialises calculator instance.
  subroutine TElectrostaticPotentials_initialise(this, input)

    !> Instance of this
    type(TElectrostaticPotentials), intent(out) :: this

    !> Input data
    type(TElectrostaticPotentialsInp), intent(inout) :: input

    this%EspOutFile = input%EspOutFile
    this%tAppendESP = input%tAppendESP
    allocate(this%ESPgrid(3,size(input%ESPgrid,dim=2)))
    this%ESPgrid = input%ESPgrid
    this%gridDimensioning = input%gridDimensioning
    this%softenESP = input%softenESP
    allocate(this%ESPpotential(size(input%ESPgrid,dim=2)))
    allocate(this%extESPpotential(size(input%ESPgrid,dim=2)))

  end subroutine TElectrostaticPotentials_initialise

  subroutine evaluate(this, env, SccCalc, EField)

    !> Object holding the potential location information
    class(TElectrostaticPotentials), intent(inout) :: this

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Module variables for SCC
    type(TScc), allocatable, intent(inout) :: sccCalc

    !> Electric field magnitude
    real(dp), intent(in) :: EField(3)

    integer :: ii

    this%ESPpotential = 0.0_dp
    this%extESPpotential = 0.0_dp

    call sccCalc%internalElectroStaticPotential(this%ESPpotential, env, this%ESPgrid,&
        & epsSoften=this%softenESP)
    call sccCalc%externalElectroStaticPotential(this%extESPpotential, env, this%ESPgrid,&
        & epsSoften=this%softenESP)

    if (any(EField /= 0.0_dp)) then
      do ii = 1, size(this%ESPgrid,dim=2)
        this%extESPpotential(ii) = this%extESPpotential(ii) + dot_product(this%ESPgrid(:, ii),&
            & EField)
      end do
    end if

  end subroutine evaluate

end module electrostaticPotentials
