!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Exporting the functionality we use from the library mbd.
module dftbp_extlibs_mbd
  use mbd, only : mbd_calc_t, TDispMbdInp => mbd_input_t
#:if WITH_MPI
  use mpi_f08, only : MPI_Comm
#:endif
  implicit none

  private
  public :: TDispMbdInp, mbd_calc_t
#:if WITH_MPI
  public :: setMbdCommunicator
#:endif


#:if WITH_MPI
  !> Sets the communicator field of the MBD input
  !!
  !! Depending on whether MBD was compiled with the mpi_f08 module or not, the communicator field is
  !! either of type(MPI_Comm) or of integer type. The generic resolves to the appropriate one.
  !!
  interface assignComm
    module procedure assignCommInt
    module procedure assignCommType
  end interface assignComm
#:endif


contains

#:if WITH_MPI

  !> Sets the MPI communicator MBD should use (otherwise MBD uses MPI_COMM_WORLD)
  subroutine setMbdCommunicator(inp, commId)

    !> MBD input
    type(TDispMbdInp), intent(inout) :: inp

    !> Fortran handle of the MPI communicator to use
    integer, intent(in) :: commId

    call assignComm(inp%comm, commId)

  end subroutine setMbdCommunicator


  !> Assigns communicator handle to an integer communicator field
  subroutine assignCommInt(comm, commId)

    !> Communicator field
    integer, intent(out) :: comm

    !> Fortran handle of the MPI communicator
    integer, intent(in) :: commId

    comm = commId

  end subroutine assignCommInt


  !> Assigns communicator handle to an mpi_f08 communicator field
  subroutine assignCommType(comm, commId)

    !> Communicator field
    type(MPI_Comm), intent(out) :: comm

    !> Fortran handle of the MPI communicator
    integer, intent(in) :: commId

    comm%mpi_val = commId

  end subroutine assignCommType

#:endif

end module dftbp_extlibs_mbd
