!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#! (TYPE, RANK, NAME) tuple for all chunk types which need to be assembled
#:set CHUNK_TYPES = [('real(dp)', 1, 'R1'), ('real(dp)', 2, 'R2'), &
    & ('real(dp)', 3, 'R3'), ('real(dp)', 4, 'R4'), &
    & ('complex(dp)', 1, 'C1'), ('complex(dp)', 2, 'C2'), &
    & ('complex(dp)', 3, 'C3'), ('complex(dp)', 4, 'C4'), &
    & ('integer', 1, 'I1')]

!> Contains routines helpful for mpi-parallelisation.
module dftbp_common_schedule
  use dftbp_common_accuracy, only : dp
  use dftbp_common_environment, only : TEnvironment
#:if WITH_MPI
  use dftbp_extlibs_mpifx, only : MPI_SUM, mpifx_allreduceip, mpifx_allgatherv
#:endif
  implicit none

  private
  public :: distributeRangeInChunks, distributeRangeInChunks2, distributeRangeWithWorkload
  public :: assembleChunks, getChunkRanges, getIndicesWithWorkload, gatherChunks
  public :: getStartAndEndIndex

#:for _, _, NAME in CHUNK_TYPES
  interface assembleChunks
    module procedure assemble${NAME}$Chunks
  end interface assembleChunks

  interface gatherChunks
    module procedure gather${NAME}$Chunks
  end interface gatherChunks  
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


  !> Distributes a range among processes within a process group
  !> and take into account that each item may have a different workload
  subroutine distributeRangeWithWorkload(env, globalFirst, globalLast, workload, indices)

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> First element of the range
    integer, intent(in) :: globalFirst

    !> Last element of the range
    integer, intent(in) :: globalLast

    !> Number of elements each item has to process
    integer, intent(in) :: workload(:)

    !> Index array to be iterated over
    integer, allocatable, intent(out) :: indices(:)

    integer :: ii

  #:if WITH_MPI
    call getIndicesWithWorkload(env%mpi%groupComm%size, env%mpi%groupComm%rank, globalFirst,&
        & globalLast, workload, indices)
  #:else
    allocate(indices(globalLast - globalFirst + 1))
    indices(:) = [(ii, ii = globalFirst, globalLast)]
  #:endif

  end subroutine distributeRangeWithWorkload


#:for DTYPE, RANK, NAME in CHUNK_TYPES

  !> Assembles the chunks by summing up contributions within a process group.
  subroutine assemble${NAME}$Chunks(env,chunks)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> Array to assemble
    ${DTYPE}$, intent(inout) :: chunks${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

  #:if WITH_MPI
    call mpifx_allreduceip(env%mpi%groupComm, chunks, MPI_SUM)
  #:endif

  end subroutine assemble${NAME}$Chunks

#:endfor

#:for DTYPE, RANK, NAME in CHUNK_TYPES

  !> Gathers chunks within a process group in a global array.
  subroutine gather${NAME}$Chunks(env, globalFirst, globalLast, chunks, composite)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> First element of the range
    integer, intent(in) :: globalFirst

    !> Last element of the range
    integer, intent(in) :: globalLast

    !> Array to gather
    ${DTYPE}$, intent(in) :: chunks${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

    !> Gathered array
    ${DTYPE}$, intent(out) :: composite${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

    integer, allocatable :: locSize(:), vOffset(:)
    integer :: localFirst, localLast, iam

  #:if WITH_MPI
    allocate(locSize(env%mpi%groupComm%size))
    allocate(vOffSet(env%mpi%groupComm%size))

    iam = env%mpi%groupComm%rank
    call getChunkRanges(env%mpi%groupComm%size, iam, globalFirst, globalLast, localFirst, localLast)
    locSize(iam + 1) = localLast - localFirst + 1
    vOffSet(iam + 1) = localFirst - 1
    
    call mpifx_allgatherv(env%mpi%globalComm, chunks, composite, locSize, vOffset)
 #:else
    composite = chunks
  #:endif

  end subroutine gather${NAME}$Chunks

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
    localLast = min(localFirst + nLocal - 1, globalLast)

  end subroutine getChunkRanges


  !> Calculate the indices for a given MPI-communicator considerung different workload
  subroutine getIndicesWithWorkload(groupSize, myRank, globalFirst, globalLast, workload, indices)

    !> Size of the group over which the chunks should be distributed
    integer, intent(in) :: groupSize

    !> Rank of the current process
    integer, intent(in) :: myRank

    !> First element of the range
    integer, intent(in) :: globalFirst

    !> Last element of the range
    integer, intent(in) :: globalLast

    !> Workload for each item
    integer, intent(in) :: workload(:)

    !> Index array to be iterated over
    integer, allocatable, intent(out) :: indices(:)

    integer :: numIndices, rank, ii
    integer, allocatable :: rankWorkload(:), indices_(:)

    allocate(indices_(globalLast - globalFirst + 1))
    allocate(rankWorkload(groupSize))

    rankWorkload(:) = 0
    indices_(:) = 0
    numIndices = 0

    do ii = globalFirst, globalLast
      rank = minloc(rankWorkload, dim=1)
      rankWorkload(rank) = rankWorkload(rank) + max(1, workload(ii))
      if (rank == myRank + 1) then
        numIndices = numIndices + 1
        indices_(numIndices) = ii
      end if
    end do

    allocate(indices(numIndices))
    indices(1:numIndices) = indices_(1:numIndices)

  end subroutine getIndicesWithWorkload


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


end module dftbp_common_schedule
