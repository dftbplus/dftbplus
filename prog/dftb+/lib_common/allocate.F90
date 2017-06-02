!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module allocation
#include "assert.h"
  implicit none

  private

  public :: ioerr
  public :: allocateError, allocateLog
  

  integer :: ioerr

  interface allocateError
    module procedure allocateError1
    module procedure allocateError2
  end interface

  integer, parameter :: nrAllocTypes = 6
  character(26) :: allocStr(nrAllocTypes) = &
      &(/'Array Allocation          ', &         ! Make sure, all strings have
      &  'Array Deallocation        ', &         ! the declared length
      &  'Pointer Allocation        ', &
      &  'Pointer Deallocation      ', &
      &  'Pointer Array Allocation  ', &
      &  'Pointer Array Deallocation' /)
  
contains  

  subroutine allocateError1(fileName, lineNr, ioerr, allocType)
    character(*), intent(in) :: fileName
    integer,      intent(in) :: lineNr
    integer,      intent(in) :: ioerr
    integer,      intent(in) :: allocType

    ASSERT((allocType >= 1) .and. (allocType <= nrAllocTypes))
    write (*, '(A,A)') "!!! FAILED :   ", trim(allocStr(allocType))
    write (*, '(A,I0)') "!!! ERR. NR.:  ", ioerr   
    write (*, '(A,A)') "!!! FILE:      ", fileName
    write (*, '(A,I0)') "!!! LINE NR.:  ", lineNr
    stop
  end subroutine allocateError1


  
  subroutine allocateError2(fileName, lineNr, msg, allocType)
    character(*), intent(in) :: fileName
    integer,      intent(in) :: lineNr
    character(*), intent(in) :: msg
    integer,      intent(in) :: allocType

    ASSERT((allocType >= 1) .and. (allocType <= nrAllocTypes))
    write (*, '(A,A)') "!!! FAILED:    ", trim(allocStr(allocType))
    write (*, '(A,A)') "!!! ERR. MSG.: ", msg
    write (*, '(A,A)') "!!! FILE:      ", fileName
    write (*, '(A,I0)') "!!! LINE NR.:  ", lineNr
    stop
  end subroutine allocateError2
  
  

  subroutine allocateLog(fileName, lineNr, allocType, allocSize)
    character(*), intent(in)           :: fileName
    integer,      intent(in)           :: lineNr
    integer,      intent(in)           :: allocType
    integer,      intent(in), optional :: allocSize
    ASSERT((allocType >= 1) .and. (allocType <= nrAllocTypes))
    write (*, '(A,A)') "Succesfull:  ", trim(allocStr(allocType))
    if (present(allocSize)) then
      write (*, '(A,I0)') "Size:        ", allocSize
    end if
    write (*, '(A,A)') "File:        ", fileName
    write (*, '(A,I0)') "Line nr.:    ", lineNr
  End Subroutine allocateLog


end module allocation
