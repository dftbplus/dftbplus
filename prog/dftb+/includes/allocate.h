!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!! -*- f90 -*-
!! Provides trapped allocate and deallocate for Fortran
!! (compile the source with DEBUG>=3 to use the verbose version)
!!

use allocation

#ifndef _allocate_h_
#define _allocate_h_

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Low-level macros.
!! Do NOT use them, use the DEBUG-level specific ones instead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define GEN_ALLOCATE_(object, objsize) \
  if (.not. allocated(object)) then;\
    allocate(object objsize, stat=ioerr);\
    if (ioerr /= 0) call allocateError(__FILE__, __LINE__, ioerr, 1);\
  else;\
    call allocateError(__FILE__, __LINE__, "Array already allocated", 1);\
  endif


#define GEN_DEALLOCATE_(object) \
  if (allocated(object)) then;\
    deallocate(object, stat=ioerr);\
    if (ioerr /= 0) call allocateError(__FILE__,__LINE__,ioerr, 2);\
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Debug turned on and its level is high enough
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if DEBUG >= 3

#define ALLOCATE_(object, objsize) \
  GEN_ALLOCATE_(object, objsize);\
  call allocateLog(__FILE__,__LINE__,1,size(object))

#define DEALLOCATE_(object) \
  GEN_DEALLOCATE_(object);\
  call allocateLog(__FILE__,__LINE__,2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Debug turned off - only check for failure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#else

#define ALLOCATE_(object, objsize) \
  GEN_ALLOCATE_(object, objsize)

#define DEALLOCATE_(object) \
  GEN_DEALLOCATE_(object)

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
