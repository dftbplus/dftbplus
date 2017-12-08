! Fortran 2008 allows to pass unallocated entities (or null pointers) for optional arguments. The
! present() query for the corresponding dummy argument should return false.
!
! Known to fail:
!  Intel 16, when checks are turned on: ifort -check  (Intel 17 works, though)
!  PGI 17.10, when checks are turned on: pgfortran -Mchkptr
!
! Workaround:
!
!   Make dummy argument allocatable instead. Note, if optional dummy argument was intent(out), it
!   must be changed to intent(inout) when converted to allocatable, otherwise it is automatically
!   deallocated at entry.
!
module mymod
  implicit none

contains

  ! Failing
  subroutine testOptional(array)
    integer, intent(in), optional :: array(:)

    if (present(array)) then
      print *, 'OPTIONAL ARRAY: ', array
    else
      print *, 'NOT PRESENT'
    end if

  end subroutine testOptional


  ! Workaround
  subroutine testAllocatable(array)
    integer, intent(in), allocatable :: array(:)

    if (allocated(array)) then
      print *, 'ALLOCATABLE ARRAY: ', array
    else
      print *, 'NOT ALLOCATED'
    end if

  end subroutine testAllocatable

end module mymod


program test
  use mymod
  implicit none

  integer, allocatable :: array(:)

  print *, 'You should see "NOT ALLOCATED" below'
  call testAllocatable(array)
  print *, 'You should see "NOT PRESENT" below'
  call testOptional(array)

end program test
