!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#! (TYPE, RANK, NAME) tuple for all chunk types which need to be assembled
#:set CHUNK_TYPES = [('real(dp)', 1, 'R1'), ('real(dp)', 2, 'R2'), ('complex(dp)', 1, 'C1'),&
    & ('complex(dp)', 2, 'C2')]

!> Contains routines helpful for mpi-parallelisation.
module schedule
#:if WITH_MPI
  use mpifx
#:endif
  use environment
  use accuracy
  implicit none
  private

  public :: distributeRangeInChunks, assembleChunks

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


end module schedule
