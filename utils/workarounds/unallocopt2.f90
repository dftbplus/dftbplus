! Intel compiler crashes if unallocated class variable is passed for optional argument
!
! Known to fail:
!   Intel 17
!
! Workaround:
!   Make dummy argument allocatable instead of optional.
!
module mymod
  implicit none

  type :: BaseType
    integer :: ii = 0
  end type BaseType

contains

  subroutine testBug(bt)
    class(BaseType), intent(inout), optional :: bt

    print *, 'PRESENT:', present(bt)
    if (present(bt)) then
      bt%ii = 42
    end if

  end subroutine testBug


  subroutine testWorkaround(bt)
    class(BaseType), intent(inout), allocatable :: bt

    print *, 'ALLOCATED:', allocated(bt)
    if (allocated(bt)) then
      bt%ii = 42
    end if

  end subroutine testWorkaround

end module mymod


program myprog
  use mymod
  implicit none

  class(BaseType), allocatable :: btClass
  type(BaseType), allocatable ::  btType

  allocate(btType)
  call move_alloc(btType, btClass)
  call testBug(btClass)
  deallocate(btClass)
  call testWorkaround(btClass)
  ! Crashes here
  call testBug(btClass)

end program myprog
