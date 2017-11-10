!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains MPI related einvironment settings
module mpienv
  use mpifx
  implicit none
  private

  public :: TMpiEnv, TMpiEnv_init


  !> Contains MPI related einvironment settings
  type :: TMpiEnv
    private

    !> Global MPI communicator
    type(mpifx_comm), public :: all

  end type TMpiEnv


contains

  !> Initializes MPI environment.
  subroutine TMpiEnv_init(this)
    type(TMpiEnv), intent(out) :: this

    ! MPI settings
    call this%all%init()

  end subroutine TMpiEnv_init


end module mpienv
