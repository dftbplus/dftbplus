!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Contains routines helpful for mpi-parallelisation.
module mpiutils
  use mpifx
  implicit none
  private

  public :: distributeRangeInChunks


contains


  !> Distributes a range in chunks over processes within a communicator
  subroutine distributeRangeInChunks(comm, globalFirst, globalLast, localFirst, localLast)

    !> MPI-communicator, over which the range should be distributed
    type(mpifx_comm), intent(in) :: comm

    !> First element of the range
    integer, intent(in) :: globalFirst

    !> Last element of the range
    integer, intent(in) :: globalLast

    !> First element to process locally
    integer, intent(out) :: localFirst

    !> Last element to process locally
    integer, intent(out) :: localLast

    integer :: rangeLength, nLocal, remainder

    rangeLength = globalLast - globalFirst + 1
    nLocal = rangeLength / comm%size
    remainder = mod(rangeLength, comm%size)
    if (comm%rank < remainder) then
      nLocal = nLocal + 1
      localFirst = globalFirst + comm%rank * nLocal
    else
      localFirst = globalFirst + remainder * (nLocal + 1) + (comm%rank - remainder) * nLocal
    end if
    localLast = min(localFirst + nLocal - 1, rangeLength)

  end subroutine distributeRangeInChunks


end module mpiutils
