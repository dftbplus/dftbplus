!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2024  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

! Repeated polymorphic assignment to an allocatable may fail, if the allocatable was transfered
! via move_alloc to another allocatable.
!
! Known to fail:
!  Intel 19, 2021
!
! https://community.intel.com/t5/Intel-Fortran-Compiler/Compiler-bug-Repeated-polymorphic-assignment-raises-run-time/td-p/1339479
!
! Workaround:
!
! Use a separate variable in each polymorphic assignment
!
program testprog
  implicit none

  type :: base_t
  end type base_t

  type, extends(base_t) :: extended_t
    integer :: ii
  end type extended_t

  class(base_t), allocatable :: base, base2, buffer

  base = extended_t(41)
  call move_alloc(base, buffer)
  !base = extended_t(42)
  base2 = extended_t(42)
  print *, "DONE"

end program testprog
