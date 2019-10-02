module xmlf90_wxml_elstack

implicit none
  
private

!
! Simple stack to keep track of which elements have appeared so far
!
integer, parameter, private            :: STACK_SIZE = 20

type, public :: elstack_t
private
      integer                                    :: n_items
      character(len=100), dimension(STACK_SIZE)  :: data
end type elstack_t

public  :: push_elstack, pop_elstack, reset_elstack, print_elstack
public  :: get_top_elstack, is_empty, get_elstack_signature
public  :: len

interface len
   module procedure number_of_items
end interface
private :: number_of_items

interface is_empty
      module procedure is_empty_elstack
end interface
private :: is_empty_elstack

CONTAINS

!-----------------------------------------------------------------
subroutine reset_elstack(elstack)
type(elstack_t), intent(inout)  :: elstack

elstack%n_items = 0

end subroutine reset_elstack

!-----------------------------------------------------------------
function is_empty_elstack(elstack) result(answer)
type(elstack_t), intent(in)  :: elstack
logical                    :: answer

answer = (elstack%n_items == 0)
end function is_empty_elstack

!-----------------------------------------------------------------
function number_of_items(elstack) result(n)
type(elstack_t), intent(in)  :: elstack
integer                      :: n

n = elstack%n_items
end function number_of_items

!-----------------------------------------------------------------
subroutine push_elstack(item,elstack)
character(len=*), intent(in)      :: item
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
character(len=*), intent(out)        :: item

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
character(len=*), intent(out)        :: item

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
      write(unit=unit,fmt=*) trim(elstack%data(i))
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
   length = len_trim(elstack%data(i))
   string(j+1:j+1) = "/"
   j = j+1
   string(j+1:j+length) = trim(elstack%data(i))
   j = j + length
enddo

end subroutine get_elstack_signature

end module xmlf90_wxml_elstack
