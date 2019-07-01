module xmlf90_dictionary

  use xmlf90_buffer

implicit none
  
private
!
! A very rough implementation for now
! It uses fixed-length buffers for key/value pairs,
! and the maximum number of dictionary items is hardwired.

integer, parameter, private    :: MAX_ITEMS = 64
type, public :: dictionary_t
private
      integer                               :: number_of_items
      type(buffer_t), dimension(MAX_ITEMS)  :: key
      type(buffer_t), dimension(MAX_ITEMS)  :: value
end type dictionary_t

!
! Building procedures
!
public  :: add_key_to_dict, add_value_to_dict, init_dict, reset_dict

!
! Query and extraction procedures
!
public  :: len
interface len
   module procedure number_of_entries
end interface
public  :: number_of_entries
public  :: get_key
public  :: get_value
public  :: has_key
public  :: print_dict
!
public  :: get_name

interface get_name
   module procedure get_key
end interface

interface get_value
   module procedure sax_get_value
end interface
private :: sax_get_value

CONTAINS

!------------------------------------------------------
function number_of_entries(dict) result(n)
type(dictionary_t), intent(in)   :: dict
integer                          :: n

n = dict%number_of_items

end function number_of_entries

!------------------------------------------------------
function has_key(dict,key) result(found)
type(dictionary_t), intent(in)   :: dict
character(len=*), intent(in)     :: key
logical                          :: found

integer  :: n, i
found = .false.
n = dict%number_of_items
do  i = 1, n
      if (dict%key(i) .EQUAL. key) then
         found = .true.
         exit
      endif
enddo
end function has_key

!------------------------------------------------------
subroutine sax_get_value(dict,key,value,status)
type(dictionary_t), intent(in)            :: dict
character(len=*), intent(in)              :: key
character(len=*), intent(out)             :: value
integer, intent(out)                      :: status
!
integer  :: n, i

status = -1
n = dict%number_of_items
do  i = 1, n
      if (dict%key(i) .EQUAL. key) then
         value = str(dict%value(i))
         status = 0
         RETURN
      endif
enddo

end subroutine sax_get_value

!------------------------------------------------------
subroutine get_key(dict,i,key,status)
!
! Get the i'th key
!
type(dictionary_t), intent(in)            :: dict
integer, intent(in)                       :: i
character(len=*), intent(out)             :: key
integer, intent(out)                      :: status

if (i <= dict%number_of_items) then
      key = str(dict%key(i))
      status = 0
else
      key = ""
      status = -1
endif

end subroutine get_key

!------------------------------------------------------
subroutine add_key_to_dict(key,dict)
type(buffer_t), intent(in)          :: key
type(dictionary_t), intent(inout)   :: dict

integer  :: n

n = dict%number_of_items
if (n == MAX_ITEMS) then
      write(unit=0,fmt=*) "Dictionary capacity exceeded ! size= ", max_items
      RETURN
endif

n = n + 1
dict%key(n) = key
dict%number_of_items = n

end subroutine add_key_to_dict

!------------------------------------------------------
! Assumes we build the dictionary in an orderly fashion,
! so one adds first the key and then immediately afterwards the value.
!
subroutine add_value_to_dict(value,dict)
type(buffer_t), intent(in)          :: value
type(dictionary_t), intent(inout)   :: dict

integer  :: n

n = dict%number_of_items
dict%value(n) = value

end subroutine add_value_to_dict

!------------------------------------------------------
subroutine init_dict(dict)
type(dictionary_t), intent(inout)   :: dict

integer  :: i

dict%number_of_items = 0
do i=1, MAX_ITEMS                      ! To avoid "undefined" status
   call init_buffer(dict%key(i))       ! (Fortran90 restriction)
   call init_buffer(dict%value(i))
enddo
end subroutine init_dict
!------------------------------------------------------
subroutine reset_dict(dict)
type(dictionary_t), intent(inout)   :: dict

dict%number_of_items = 0

end subroutine reset_dict

!------------------------------------------------------
subroutine print_dict(dict)
type(dictionary_t), intent(in)   :: dict

integer  :: i

do i = 1, dict%number_of_items
      print *, trim(str(dict%key(i))), " = ", trim(str(dict%value(i)))
enddo

end subroutine print_dict


end module xmlf90_dictionary
