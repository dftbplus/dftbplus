module xmlf90_elstack

use xmlf90_buffer

implicit none

private

!
! Simple stack to keep track of which elements have appeared so far
!
integer, parameter, private            :: STACK_SIZE = 40

type, public :: elstack_t
private
      integer                                :: n_items
      type(buffer_t), dimension(STACK_SIZE)  :: data
end type elstack_t

public  :: push_elstack, pop_elstack, reset_elstack, print_elstack
public  :: init_elstack
public  :: get_top_elstack, is_empty, get_elstack_signature

interface is_empty
      module procedure is_empty_elstack
end interface
private :: is_empty_elstack

CONTAINS

!-----------------------------------------------------------------
subroutine init_elstack(elstack)
type(elstack_t), intent(inout)  :: elstack

integer :: i

elstack%n_items = 0
do i = 1, STACK_SIZE                   ! to avoid "undefined status"
      call init_buffer(elstack%data(i))
enddo
end subroutine init_elstack

!-----------------------------------------------------------------
subroutine reset_elstack(elstack)
type(elstack_t), intent(inout)  :: elstack

integer :: i

elstack%n_items = 0
do i = 1, STACK_SIZE                  
      call reset_buffer(elstack%data(i))
enddo
end subroutine reset_elstack

!-----------------------------------------------------------------
function is_empty_elstack(elstack) result(answer)
type(elstack_t), intent(in)  :: elstack
logical                    :: answer

answer = (elstack%n_items == 0)
end function is_empty_elstack

!-----------------------------------------------------------------
subroutine push_elstack(item,elstack)
type(buffer_t), intent(in)      :: item
type(elstack_t), intent(inout)  :: elstack

integer   :: n

n = elstack%n_items
if (n == STACK_SIZE) then
      stop "*Element stack full"
endif
n = n + 1
elstack%data(n) = item
elstack%n_items = n

end subroutine push_elstack

!-----------------------------------------------------------------
subroutine pop_elstack(elstack,item)
type(elstack_t), intent(inout)     :: elstack
type(buffer_t), intent(out)        :: item

!
! We assume the elstack is not empty... (the user has called is_empty first)
!
integer   :: n

n = elstack%n_items
if (n == 0) then
      stop "*********Element stack empty"
endif
item = elstack%data(n)
elstack%n_items = n - 1

end subroutine pop_elstack

!-----------------------------------------------------------------
subroutine get_top_elstack(elstack,item)
!
! Get the top element of the stack, *without popping it*.
!
type(elstack_t), intent(in)        :: elstack
type(buffer_t), intent(out)        :: item

!
! We assume the elstack is not empty... (the user has called is_empty first)
!
integer   :: n

n = elstack%n_items
if (n == 0) then
      stop "*********Element stack empty"
endif
item = elstack%data(n)

end subroutine get_top_elstack

!-----------------------------------------------------------------
subroutine print_elstack(elstack,unit)
type(elstack_t), intent(in)   :: elstack
integer, intent(in)           :: unit
integer   :: i

do i = elstack%n_items, 1, -1
      write(unit=unit,fmt=*) str(elstack%data(i))
enddo

end subroutine print_elstack

!-------------------------------------------------------------
subroutine get_elstack_signature(elstack,string)
type(elstack_t), intent(in)   :: elstack
character(len=*), intent(out) :: string
integer   :: i, length, j

string = ""
j = 0
do i = 1, elstack%n_items
   length = len(elstack%data(i))
   string(j+1:j+1) = "/"
   j = j+1
   string(j+1:j+length) = str(elstack%data(i))
   j = j + length
enddo

end subroutine get_elstack_signature

end module xmlf90_elstack





