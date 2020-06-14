!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> IDs for file operations
module dftbp_fileid
  use dftbp_message

  private

  public :: getFileId


  !> starting range of IDs
  integer, parameter :: minId = 20

  !> largest ID
  integer, parameter :: maxId = 65535

contains


  !> get a new (unused) file ID
  function getFileId()
    integer :: getFileId

    integer, save :: curId = minId

    if (curId > maxId) then
      call error("getFileId: No more free file identification numbers")
    end if
    getFileId = curId
    curId = curId + 1

  end function getFileId

end module dftbp_fileid
