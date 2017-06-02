!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!! -*- f90 -*-
!! vim:syntax=fortran:
!!
!! Provides C-style assertions for Fortran
!! (compile the source with DEBUG>=1 to use it)
!!

#if DEBUG > 0
use assert
#endif

#ifndef _assert_h_
#define _assert_h_


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Debug turned on
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if DEBUG > 0

#define ASSERT(condition) \
  if (.not.(condition)) call assertError(__FILE__, __LINE__)

#define EMB_ASSERT(condition, file, line) \
  if (.not.(condition)) call assertError(file, line)

#define ASSERT_ENV(statement) \
  statement


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Debug turned off
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#else

#define ASSERT(condition)

#define EMB_ASSERT(condition, file, line)

#define ASSERT_ENV(statement)

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
