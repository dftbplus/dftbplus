!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides auxiliary routines for parallelization.
module dftbp_common_parallel
  use dftbp_common_environment, only : TEnvironment

  implicit none
  private

  public :: getStartAndEndIndex


contains

  !> Returns the start and end index of an MPI process that calculates parts of a loop.
  subroutine getStartAndEndIndex(env, nElements, iStart, iEnd)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Array size to split
    integer, intent(in) :: nElements

    !> Start and end index of current element range
    integer, intent(out) :: iStart, iEnd

  #:if WITH_MPI
    !! Size of split index regions
    integer :: splitSize

    !! Number of elements that exceed integer times nProcs
    integer :: offset
  #:endif

    @:ASSERT(nElements >= 0)

  #:if WITH_MPI
    @:ASSERT(env%mpi%globalComm%rank < env%mpi%globalComm%size)

    splitSize = nElements / env%mpi%globalComm%size

    ! start and end indices assuming equal split sizes
    iStart = env%mpi%globalComm%rank * splitSize + 1
    iEnd = iStart + splitSize - 1

    ! distribute possible remainder to the ranges at the end
    offset = env%mpi%globalComm%size - mod(nElements, env%mpi%globalComm%size)
    if (env%mpi%globalComm%rank + 1 > offset) then
      iStart = iStart + env%mpi%globalComm%rank - offset
      iEnd = iEnd + env%mpi%globalComm%rank - offset + 1
    end if
  #:else
    iStart = 1
    iEnd = nElements
  #:endif

  end subroutine getStartAndEndIndex

end module dftbp_common_parallel
