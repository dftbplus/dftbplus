!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************

#:include "error.fypp"

#!#:set MEMLOG = 1

#! (LABEL, TYPE, ARRAY, ARGS) tuple for all logged arrays
#:set ALLOC_CASES = [&
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
  use, intrinsic :: iso_fortran_env, only : int64
  use, intrinsic :: iso_c_binding, only : c_sizeof
  use std_io
  use dftbp_accuracy, only : dp

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
  subroutine allocate_${LABEL}$(array, ${ARGS}$, err)

    ${TYPE}$, allocatable, target, intent(inout) :: array(${ARRAY}$)

    integer, intent(in) :: ${ARGS}$

    !> Error code, 0 if no problems
    integer, intent(out), optional :: err

    ${TYPE}$, pointer :: pArrayFlat(:)
    integer :: iErr

    if (present(err)) then
      err = 0
    end if

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
      @:ERROR_HANDLING(err,-1,"ALLOCATION ERROR: ${TYPE}$ (${LABEL}$) array is already allocated")
    endif

    if(.not. allocated(array)) then
      allocate(array(${ARGS}$), stat=iErr)
      if (ierr /= 0) then
        @:FORMATTED_ERROR_HANDLING(err, iErr, "(A,I0)", "Poisson allocation error: ", iErr)
      else
        if (size(array) > 0) then
          pArrayFlat(1 : size(array)) => array
          alloc_mem= alloc_mem + c_sizeof(pArrayFlat(1)) * size(pArrayFlat)
          if (alloc_mem > peak_mem) then
            peak_mem = alloc_mem
          endif
        end if
      #:if defined('MEMLOG')
        call writePoissMemInfo()
      #:endif
      endif
    endif

  end subroutine allocate_${LABEL}$

#:endfor


#:for LABEL, TYPE, ARRAY, _ in ALLOC_CASES

  !---------------------------------------------------------------
  subroutine deallocate_${LABEL}$(array, err)

    ${TYPE}$, allocatable, target, intent(inout) :: array(${ARRAY}$)

    !> Error code, 0 if no problems
    integer, intent(out), optional :: err

    ${TYPE}$, pointer :: pArrayFlat(:)

    if (present(err)) then
      err = 0
    end if

    if (allocated(array)) then
      if (size(array) > 0) then
        pArrayFlat(1 : size(array)) => array
        alloc_mem= alloc_mem - c_sizeof(pArrayFlat(1)) * size(pArrayFlat)
      end if
      deallocate(array)
     #:if defined('MEMLOG')
      call writePoissMemInfo()
     #:endif
    else
      @:ERROR_HANDLING(err, -1, 'ALLOCATION ERROR: ${TYPE}$ (${LABEL}$) array is not allocated')
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

    if (mem < 1000_int64) then
      str = ' bt'
      dec = 1.0_dp
    else if (mem < 10000000_int64) then
      str = ' kb'
      dec = 1.0E-3_dp
    else if (mem < 10000000000_int64) then
      str = ' Mb'
      dec = 1.0E-6_dp
    else
      str = ' Gb'
      dec = 1.0E-9_dp
    endif

  end subroutine memstr

end module gallocation
