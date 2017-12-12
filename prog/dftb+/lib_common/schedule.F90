!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

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

#:for NAME in [('R1'), ('R2'), ('C1'), ('C2')]
  interface assembleChunks
    module procedure assemble${NAME}$Chunks
  end interface assembleChunks
#:endfor

contains


  !> Distributes a range in chunks over processes within a communicator
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

    integer :: rangeLength, nLocal, remainder

    rangeLength = globalLast - globalFirst + 1
    nLocal = rangeLength / env%mpi%groupComm%size
    remainder = mod(rangeLength, env%mpi%groupComm%size)
    if (env%mpi%groupComm%rank < remainder) then
      nLocal = nLocal + 1
      localFirst = globalFirst + env%mpi%groupComm%rank * nLocal
    else
      localFirst = globalFirst + remainder * (nLocal + 1)&
          & + (env%mpi%groupComm%rank - remainder) * nLocal
    end if
    localLast = min(localFirst + nLocal - 1, rangeLength)

#:else

    localFirst = globalFirst
    localLast = globalLast

#:endif

  end subroutine distributeRangeInChunks


#:for DTYPE, NAME, SHAPE in [('real(dp)', 'R1', '(:)'), ('real(dp)', 'R2', '(:,:)'), ('complex(dp)', 'C1', '(:)'), ('complex(dp)', 'C2', '(:,:)')]

  subroutine assemble${NAME}$Chunks(env,chunks)

    !> Environment settings
    type(TEnvironment), intent(in) :: env

    !> array to assemble
    ${DTYPE}$, intent(inout) :: chunks${SHAPE}$

#:if WITH_MPI
    call mpifx_allreduceip(env%mpi%groupComm, chunks, MPI_SUM)
#:endif

  end subroutine assemble${NAME}$Chunks

#:endfor

end module schedule
