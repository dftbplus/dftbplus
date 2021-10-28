!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module waveplot_parallel

  use dftbp_common_accuracy, only : dp

  implicit none

  private

  public :: getStartAndEndIndex, getStartAndEndIndices


contains

  pure subroutine getStartAndEndIndex(nSystems, nProcs, iProc, iStart, iEnd)

    !> array size to split
    integer, intent(in) :: nSystems

    !> number of available processes
    integer, intent(in) :: nProcs

    !> current process index
    integer, intent(in) :: iProc

    !> start and end index of current chunk
    integer, intent(out) :: iStart, iEnd

    !> size of splitted index regions
    integer :: splitSize

    !> number of systems that exceeds integer times nProcs
    integer :: offset

    splitSize = nSystems / nProcs

    ! start and end indices assuming equal split sizes
    iStart = iProc * splitSize + 1
    iEnd = iStart + splitSize - 1

    ! distribute possible remainder to the chunks at the end
    offset = nProcs - mod(nSystems, nProcs)
    if (iProc + 1 > offset) then
      iStart = iStart + iProc - offset
      iEnd = iEnd + iProc - offset + 1
    end if

  end subroutine getStartAndEndIndex


  pure subroutine getStartAndEndIndices(nSystems, nChunks, chunks)

    !> array size to split
    integer, intent(in) :: nSystems

    !> number of chunks to build
    integer, intent(in) :: nChunks

    !> start and end indices of all chunks
    integer, intent(out), allocatable :: chunks(:,:)

    !> auxiliary variable
    integer :: iChunk

    allocate(chunks(2, nChunks))

    do iChunk = 1, nChunks
      call getStartAndEndIndex(nSystems, nChunks, iChunk - 1, chunks(1, iChunk), chunks(2, iChunk))
    end do

  end subroutine getStartAndEndIndices

end module waveplot_parallel
