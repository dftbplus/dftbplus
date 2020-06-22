!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#! (TYPE, RANK, NAME) tuple for all chunk types which need to be assembled
#:set CHUNK_TYPES = [('real(dp)', 1, 'R1'), ('real(dp)', 2, 'R2'), ('real(dp)', 3, 'R3'), &
    & ('complex(dp)', 1, 'C1'), ('complex(dp)', 2, 'C2'), ('complex(dp)', 3, 'C3')]

!> Contains routines helpful for mpi-parallelisation.
module dftbp_schedule
#:if WITH_MPI
  use dftbp_mpifx
#:endif
  use dftbp_environment
  use dftbp_accuracy
  implicit none
  private

  public :: distributeRangeInChunks, distributeRangeInChunks2
  public :: assembleChunks

#:for _, _, NAME in CHUNK_TYPES
  interface assembleChunks
    module procedure assemble${NAME}$Chunks
  end interface assembleChunks
#:endfor

contains


  !> Distributes a range in chunks over processes within a process group.
  subroutine distributeRangeInChunks(env, globalFirst, globalLast, localFirst, localLast)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> First element of the range
    integer, intent(in) :: globalFirst

    !> Last element of the range
    integer, intent(in) :: globalLast

    !> First element to process locally
    integer, intent(out) :: localFirst

    !> Last element to process locally
    integer, intent(out) :: localLast

  #:if WITH_MPI
    call getChunkRanges(env%mpi%groupComm%size, env%mpi%groupComm%rank, globalFirst, globalLast,&
        & localFirst, localLast)
  #:else
    localFirst = globalFirst
    localLast = globalLast
  #:endif

  end subroutine distributeRangeInChunks


  !> Distributes a ranges in a double loop in chunks over processes within a process group.
  !>
  !> It will chop the loop with the wider range into chunks and leave the other intact.
  !>
  subroutine distributeRangeInChunks2(env, globalFirst1, globalLast1, globalFirst2, globalLast2,&
      & localFirst1, localLast1, localFirst2, localLast2)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> First element of the range for the outer loop
    integer, intent(in) :: globalFirst1

    !> Last element of the range for the outer loop
    integer, intent(in) :: globalLast1

    !> First element of the range for the inner loop
    integer, intent(in) :: globalFirst2

    !> Last element of the range for the inner loop
    integer, intent(in) :: globalLast2

    !> First element to process locally
    integer, intent(out) :: localFirst1

    !> Last element to process locally
    integer, intent(out) :: localLast1

    !> First element to process locally
    integer, intent(out) :: localFirst2

    !> Last element to process locally
    integer, intent(out) :: localLast2

  #:if WITH_MPI
    if (globalLast1 - globalFirst1 >= globalLast2 - globalFirst2) then
      call getChunkRanges(env%mpi%groupComm%size, env%mpi%groupComm%rank, globalFirst1,&
          & globalLast1, localFirst1, localLast1)
      localFirst2 = globalFirst2
      localLast2 = globalLast2
    else
      localFirst1 = globalFirst1
      localLast1 = globalLast1
      call getChunkRanges(env%mpi%groupComm%size, env%mpi%groupComm%rank, globalFirst2,&
          & globalLast2, localFirst2, localLast2)
    end if
  #:else
    localFirst1 = globalFirst1
    localLast1 = globalLast1
    localFirst2 = globalFirst2
    localLast2 = globalLast2
  #:endif

  end subroutine distributeRangeInChunks2


#:for DTYPE, RANK, NAME in CHUNK_TYPES

  !> Assembles the chunks by summing up contributions within a process group.
  subroutine assemble${NAME}$Chunks(env,chunks)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> array to assemble
    ${DTYPE}$, intent(inout) :: chunks${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

  #:if WITH_MPI
    call mpifx_allreduceip(env%mpi%groupComm, chunks, MPI_SUM)
  #:endif

  end subroutine assemble${NAME}$Chunks

#:endfor


  !> Calculate the chunk ranges for a given MPI-communicator.
  subroutine getChunkRanges(groupSize, myRank, globalFirst, globalLast, localFirst, localLast)

    !> Size of the group over which the chunks should be distributed
    integer, intent(in) :: groupSize

    !> Rank of the current process
    integer, intent(in) :: myRank

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
    nLocal = rangeLength / groupSize
    remainder = mod(rangeLength, groupSize)
    if (myRank < remainder) then
      nLocal = nLocal + 1
      localFirst = globalFirst + myRank * nLocal
    else
      localFirst = globalFirst + remainder * (nLocal + 1) + (myRank - remainder) * nLocal
    end if
    localLast = min(localFirst + nLocal - 1, rangeLength)

  end subroutine getChunkRanges


end module dftbp_schedule
