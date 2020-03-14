!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************

#!#:set MEMLOG = 1

#! (LABEL, TYPE, ARRAY, ARGS) tuple for all logged arrays
#:set ALLOC_CASES = [('l','logical',':', 'length'),&
  & ('i', 'integer', ':', 'length'),&
  & ('d', 'real(dp)', ':', 'length'),&
  & ('z', 'complex(dp)', ':', 'length'),&
  & ('i2', 'integer', ':,:', 'row,col'),&
  & ('d2', 'real(dp)', ':,:', 'row,col'),&
  & ('z2', 'complex(dp)', ':,:', 'row,col'),&
  & ('i3', 'integer', ':,:,:', 'row,col,dep'),&
  & ('d3', 'real(dp)', ':,:,:', 'row,col,dep'),&
  & ('z3', 'complex(dp)', ':,:,:', 'row,col,dep'),&
  & ('d4', 'real(dp)', ':,:,:,:', 'row,col,dep,qep')]

module gallocation
  use std_io
  use dftbp_accuracy, only : dp, lc
  use dftbp_message
  use, intrinsic :: iso_fortran_env, only : int64

  integer, parameter :: long = int64
  integer, save :: iolog
  integer(int64), save :: alloc_mem, peak_mem

  public log_gallocate, log_gdeallocate, writePoissMemInfo, writePoissPeakInfo

  interface log_gallocate
  #:for LABEL, _, _, _ in ALLOC_CASES
    module procedure allocate_${LABEL}$
  #:endfor
  end interface

  interface log_gdeallocate
  #:for LABEL, _, _, _ in ALLOC_CASES
    module procedure deallocate_${LABEL}$
  #:endfor
  end interface

  !---------------------------------------------------------------
contains

#:for LABEL, TYPE, ARRAY, ARGS in ALLOC_CASES

  !---------------------------------------------------------------
  subroutine allocate_${LABEL}$(array, ${ARGS}$)

    ${TYPE}$, allocatable, intent(inout) :: array(${ARRAY}$)

    integer, intent(in) :: ${ARGS}$

    integer :: iErr
    character(lc) :: strTmp

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
      call error('ALLOCATION ERROR: ${TYPE}$ (${LABEL}$) array is already allocated')
    endif

    if(.not. allocated(array)) then
      allocate(array(${ARGS}$), stat=iErr)
      if (ierr /= 0) then
        write(strTmp, "(A,I0)")'Poisson allocation error: ', iErr
        call error(strTmp)
      else
        alloc_mem= alloc_mem + sizeof(array)
        if (alloc_mem > peak_mem) then
          peak_mem = alloc_mem
        endif
      #:if defined('MEMLOG')
        call writePoissMemInfo()
      #:endif
      endif
    endif

  end subroutine allocate_${LABEL}$

#:endfor


#:for LABEL, TYPE, ARRAY, _ in ALLOC_CASES

  !---------------------------------------------------------------
  subroutine deallocate_${LABEL}$(array)

    ${TYPE}$, allocatable, intent(inout) :: array(${ARRAY}$)

    if (allocated(array)) then
      alloc_mem= alloc_mem - sizeof(array)
      deallocate(array)
     #:if defined('MEMLOG')
      call writePoissMemInfo()
     #:endif
    else
      call error('ALLOCATION ERROR: ${TYPE}$ (${LABEL}$) array is not allocated')
    endif

  end subroutine deallocate_${LABEL}$

#:endfor


  ! ------------------------------------------------------------
  subroutine writePoissMemInfo(fileId)

    integer, intent(in), optional :: fileId

    character(3) :: str
    real(dp) :: dec

    call memstr(alloc_mem,dec,str)
    if (present(fileId)) then
      write(fileId,'(A35,F8.2,A3)') 'current Poisson memory allocated: ', real(alloc_mem,dp)*dec,&
          & str
    else
      write(stdOut,'(A35,F8.2,A3)') 'current Poisson memory allocated: ', real(alloc_mem,dp)*dec,&
          & str
    end if

  end subroutine writePoissMemInfo

  ! ------------------------------------------------------------
  subroutine writePoissPeakInfo(fileId)

    integer, intent(in), optional :: fileId

    character(3) :: str
    real(dp) :: dec

    call memstr(peak_mem,dec,str)
    if (present(fileId)) then
      write(fileId,'(A32,T36,F8.2,A3)') 'peak Poisson memory allocated: ', real(peak_mem,dp)*dec,&
          & str
    else
      write(stdOut,'(A32,T36,F8.2,A3)') 'peak Poisson memory allocated: ', real(peak_mem,dp)*dec,&
          & str
    end if

  end subroutine writePoissPeakInfo


  !> Labels for different sizes of value
  pure subroutine memstr(mem,dec,str)

    integer(int64), intent(in) :: mem
    real(dp), intent(out) :: dec
    character(3), intent(out) :: str

    if (mem < 1000) then
      str = ' bt'
      dec = 1.0_dp
    else if (mem < 1E7_int64) then
      str = ' kb'
      dec = 1.0E-3_dp
    else if (mem < 1E10_int64) then
      str = ' Mb'
      dec = 1.0E-6_dp
    else
      str = ' Gb'
      dec = 1.0E-9_dp
    endif

  end subroutine memstr

end module gallocation
