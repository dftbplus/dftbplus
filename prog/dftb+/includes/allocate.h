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


#define GEN_INIT_P(ptr) \
  nullify(ptr)


#define GEN_ALLOCATE_P(ptr) \
  if (.not. associated(ptr)) then;\
    allocate(ptr, stat=ioerr);\
    if (ioerr /= 0) call allocateError(__FILE__, __LINE__, ioerr, 3);\
  else;\
    call allocateError(__FILE__, __LINE__, "Ptr. already associated!", 3);\
  endif


#define GEN_DEALLOCATE_P(ptr) \
  if (associated(ptr)) then;\
    deallocate(ptr, stat=ioerr);\
    if (ioerr /= 0) call allocateError(__FILE__, __LINE__, ioerr, 4);\
    GEN_INIT_P(ptr);\
  endif


#define GEN_INIT_PARR(ptr) \
  nullify(ptr)


#define GEN_ALLOCATE_PARR(ptr, dim) \
  if (.not. associated(ptr)) then;\
    allocate(ptr dim, stat=ioerr);\
    if (ioerr /= 0) call allocateError(__FILE__, __LINE__, ioerr, 5);\
  else;\
    call allocateError(__FILE__, __LINE__, "Ptr. array already associated!",5);\
  endif


#define GEN_DEALLOCATE_PARR(ptr) \
  if (associated(ptr)) then;\
    deallocate(ptr, stat=ioerr);\
    if (ioerr /= 0) call allocateError(__FILE__, __LINE__, ioerr, 6);\
    GEN_INIT_PARR(ptr);\
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

#define INIT_P(ptr) \
  GEN_INIT_P(ptr)

#define ALLOCATE_P(ptr) \
  GEN_ALLOCATE_P(ptr);\
  call allocateLog(__FILE__, __LINE__,3)

#define DEALLOCATE_P(ptr) \
  GEN_DEALLOCATE_P(ptr);\
  call allocateLog(__FILE__, __LINE__,4)

#define INIT_PARR(ptr) \
  GEN_INIT_PARR(ptr)

#define ALLOCATE_PARR(ptr, dim) \
  GEN_ALLOCATE_PARR(ptr, dim);\
  call allocateLog(__FILE__, __LINE__,5,size(ptr))

#define DEALLOCATE_PARR(ptr) \
  GEN_DEALLOCATE_PARR(ptr);\
  call allocateLog(__FILE__, __LINE__,6)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Debug turned off - only check for failure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#else

#define ALLOCATE_(object, objsize) \
  GEN_ALLOCATE_(object, objsize)

#define DEALLOCATE_(object) \
  GEN_DEALLOCATE_(object)

#define INIT_P(ptr) \
  GEN_INIT_P(ptr)

#define ALLOCATE_P(ptr) \
  GEN_ALLOCATE_P(ptr)

#define DEALLOCATE_P(ptr) \
  GEN_DEALLOCATE_P(ptr)

#define INIT_PARR(ptr) \
  GEN_INIT_PARR(ptr)

#define ALLOCATE_PARR(ptr, dim) \
  GEN_ALLOCATE_PARR(ptr, dim)

#define DEALLOCATE_PARR(ptr) \
  GEN_DEALLOCATE_PARR(ptr)

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Macros for initialization and allocation in one step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define INITALLOCATE_P(ptr) \
  INIT_P(ptr);\
  ALLOCATE_P(ptr)

#define INITALLOCATE_PARR(ptr, dim) \
  INIT_PARR(ptr);\
  ALLOCATE_PARR(ptr, dim)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
