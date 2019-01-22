!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************
!!$#define "MEMLOG"

module gallocation
  use std_io
  use gprecision

  integer, save :: iolog
  INTEGER(long), SAVE :: alloc_mem, peak_mem

  public log_gallocate, log_gdeallocate, writeMemInfo, writePeakInfo
  public resetMemLog, openMemLog, writeMemLog, closeMemLog

  interface log_gallocatep
    module procedure allocate_pd, allocate_pi, allocate_pz
    module procedure allocate_pd2, allocate_pi2, allocate_pz2
  end interface

  interface log_gallocate
     module procedure allocate_l
     module procedure allocate_d, allocate_i, allocate_z
     module procedure allocate_d2, allocate_i2, allocate_z2
     module procedure allocate_d3, allocate_i3, allocate_z3
     module procedure allocate_d4
  end interface

  interface log_gdeallocatep
    module procedure deallocate_pd, deallocate_pi, deallocate_pz
    module procedure deallocate_pd2, deallocate_pi2, deallocate_pz2
  end interface

  interface log_gdeallocate
     module procedure deallocate_l
     module procedure deallocate_d, deallocate_i, deallocate_z
     module procedure deallocate_d2, deallocate_i2, deallocate_z2
     module procedure deallocate_d3, deallocate_i3, deallocate_z3
     module procedure deallocate_d4
  end interface


  !---------------------------------------------------------------
  !---------------------------------------------------------------
contains
  !---------------------------------------------------------------
  subroutine allocate_pi(array,length)
    integer, DIMENSION(:), POINTER :: array
    integer :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_pi

  subroutine allocate_pd(array,length)
    real(kind=dp), DIMENSION(:), POINTER :: array
    integer :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_pd

  subroutine allocate_pz(array,length)
    complex(kind=dp), DIMENSION(:), POINTER :: array
    integer :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*2*dp
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_pz

    !---------------------------------------------------------------
  subroutine allocate_pi2(array,row,col)
    integer, DIMENSION(:,:), POINTER :: array
    integer :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_pi2

  subroutine allocate_pd2(array,row,col)
    real(kind=dp), DIMENSION(:,:), POINTER :: array
    integer :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_pd2

  subroutine allocate_pz2(array,row,col)
    complex(kind=dp), DIMENSION(:,:), POINTER :: array
    integer :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (associated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. associated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*2*dp
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_pz2
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  subroutine allocate_i(array,length)
    integer, DIMENSION(:), ALLOCATABLE :: array
    integer :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_i

  subroutine allocate_d(array,length)
    real(kind=dp), DIMENSION(:), ALLOCATABLE :: array
    integer :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (ALLOCATED(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not.ALLOCATED(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_d

  subroutine allocate_z(array,length)
    complex(kind=dp), DIMENSION(:), ALLOCATABLE :: array
    integer :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (ALLOCATED(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not.ALLOCATED(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*2*dp
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_z
  !---------------------------------------------------------------

  subroutine allocate_i2(array,row,col)
    integer, DIMENSION(:,:), ALLOCATABLE :: array
    integer :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_i2

  subroutine allocate_d2(array,row,col)
    real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array
    integer :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_d2

  subroutine allocate_z2(array,row,col)
    complex(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array
    integer :: row,col,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*2*dp
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_z2
  !---------------------------------------------------------------

  subroutine allocate_i3(array,row,col,dep)
    integer, DIMENSION(:,:,:), ALLOCATABLE :: array
    integer :: row,col,dep,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col,dep),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_i3

  subroutine allocate_d3(array,row,col,dep)
    real(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array
    integer :: row,col,dep,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col,dep),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_d3

  subroutine allocate_z3(array,row,col,dep)
    complex(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array
    integer :: row,col,dep,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col,dep),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*2*dp
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_z3

  subroutine allocate_d4(array,row,col,dep,qep)
    real(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: array
    integer :: row,col,dep,qep,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(row,col,dep,qep),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*dp
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_d4

  !---------------------------------------------------------------
  !---------------------------------------------------------------
  subroutine allocate_l(array,length)
    logical, DIMENSION(:), ALLOCATABLE :: array
    integer :: length,ierr

    !Allocation control: if array is already allocated STOP and write error statement
    if (allocated(array)) then
       STOP 'ALLOCATION ERROR: array is already allocated'
    endif

    if(.not. allocated(array)) then
       allocate(array(length),stat=ierr)
       if (ierr.ne.0) then
          write(stdOut,*) "ALLOCATION ERROR"; STOP
       else
          alloc_mem= alloc_mem + size(array)*4
          if (alloc_mem.gt.peak_mem) then
             peak_mem = alloc_mem
          endif
#:if defined("MEMLOG")
          call writeMemLog
#:endif
       endif
    endif
  end subroutine allocate_l

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  subroutine deallocate_pi(array)
    integer, DIMENSION(:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_pi

  subroutine deallocate_pd(array)
    real(kind=dp), DIMENSION(:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*dp
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif

    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_pd

  subroutine deallocate_pz(array)
    complex(kind=dp), DIMENSION(:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*2*dp
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_pz
  !---------------------------------------------------------------

  subroutine deallocate_pi2(array)
    integer, DIMENSION(:,:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_pi2

  subroutine deallocate_pd2(array)
    real(kind=dp), DIMENSION(:,:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*dp
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif

    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_pd2

  subroutine deallocate_pz2(array)
    complex(kind=dp), DIMENSION(:,:), POINTER :: array

    if (associated(array)) then
       alloc_mem= alloc_mem - size(array)*2*dp
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
     endif
   end subroutine deallocate_pz2


  !---------------------------------------------------------------
  !---------------------------------------------------------------
  subroutine deallocate_l(array)
    logical, DIMENSION(:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_l

  subroutine deallocate_i(array)
    integer, DIMENSION(:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_i

  subroutine deallocate_d(array)
    real(kind=dp), DIMENSION(:), allocatable :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*dp
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif

    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_d

  subroutine deallocate_z(array)
    complex(kind=dp), DIMENSION(:), allocatable :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*2*dp
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_z

  ! ------------------------------------------------------------
  subroutine deallocate_i2(array)
    integer, DIMENSION(:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_i2

  subroutine deallocate_d2(array)
    real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*dp
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_d2

  subroutine deallocate_z2(array)
    complex(kind=dp), DIMENSION(:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*2*dp
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_z2
  ! ------------------------------------------------------------

  subroutine deallocate_i3(array)
    integer, DIMENSION(:,:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*4
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_i3

  subroutine deallocate_d3(array)
    real(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*dp
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_d3

  subroutine deallocate_z3(array)
    complex(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*2*dp
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_z3

  subroutine deallocate_d4(array)
    real(kind=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: array

    if (allocated(array)) then
       alloc_mem= alloc_mem - size(array)*dp
       deallocate(array)
#:if defined("MEMLOG")
       call writeMemLog
#:endif
    else
       write(stdOut,*) 'Warning in deallocation: array is not allocated'
    endif
  end subroutine deallocate_d4

  ! ------------------------------------------------------------
  subroutine resetMemLog
     alloc_mem=0
     peak_mem=0
  end subroutine resetMemLog
  ! ------------------------------------------------------------
  subroutine writeMemInfo(iofile)

    integer iofile
    character(3) :: str
    integer :: dec

    call memstr(alloc_mem,dec,str)
    write(iofile,'(A26,F8.2,A3)') 'current memory allocated: ',real(alloc_mem)*1.0/dec,str

  end subroutine writeMemInfo
  ! ------------------------------------------------------------
  subroutine writePeakInfo(iofile)

    integer iofile
    character(3) :: str
    integer :: dec

    call memstr(peak_mem,dec,str)
    write(iofile,'(A26,F8.2,A3)') 'peak memory allocated: ',real(peak_mem)*1.0/dec,str

  end subroutine writePeakInfo

  ! ------------------------------------------------------------
  subroutine openMemLog(iofile)
    integer iofile,err

    if(iofile.ne.6) then
       open(iofile,file='memory.log',iostat=err)
       if (err.ne.0) then
          write(stdOut,*) 'Cannot open memory log-file'
          stop
       endif
    endif
    iolog=iofile

  end subroutine openMemLog

  ! ------------------------------------------------------------
  subroutine writeMemLog
    call writeMemInfo(iolog)
  end subroutine writeMemLog

  ! ------------------------------------------------------------
  subroutine closeMemLog
    call writePeakInfo(iolog)
    close(iolog)
  end subroutine closeMemLog

  ! ------------------------------------------------------------
  subroutine memstr(mem,dec,str)

  character(3) :: str
  integer(long) :: mem
integer :: dec

    if(mem.lt.1000) then
       str=' bt'; dec=1
       return
    endif

    if(mem.lt.10000000) then
       str=' kb'; dec=1000
       return
    endif

    if(mem.ge.10000000) then
       str=' Mb'; dec=1000000
       return
    endif

  end subroutine memstr

end module gallocation
