module xmlf90_wxml_dictionary

implicit none
  
private
!
! A very rough implementation for now
! It uses fixed-length buffers for key/value pairs,
! and the maximum number of dictionary items is hardwired.

integer, parameter, private    :: MAX_ITEMS = 30
type, public :: wxml_dictionary_t
private
      integer                               :: number_of_items ! = 0
      character(len=100), dimension(MAX_ITEMS)  :: key
      character(len=100), dimension(MAX_ITEMS)  :: value
end type wxml_dictionary_t

!
! Building procedures
!
public  :: add_key_to_dict, add_value_to_dict, reset_dict

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
interface get_value
   module procedure wxml_get_value
end interface

CONTAINS

!------------------------------------------------------
function number_of_entries(dict) result(n)
type(wxml_dictionary_t), intent(in)   :: dict
integer                          :: n

n = dict%number_of_items

end function number_of_entries

!------------------------------------------------------
function has_key(dict,key) result(found)
type(wxml_dictionary_t), intent(in)   :: dict
character(len=*), intent(in)     :: key
logical                          :: found

integer  :: n, i
found = .false.
n = dict%number_of_items
do  i = 1, n
      if (dict%key(i) == key) then
         found = .true.
         exit
      endif
enddo
end function has_key

!------------------------------------------------------
subroutine wxml_get_value(dict,key,value,status)
type(wxml_dictionary_t), intent(in)            :: dict
character(len=*), intent(in)              :: key
character(len=*), intent(out)             :: value
integer, intent(out)                      :: status
!
integer  :: n, i

status = -1
n = dict%number_of_items
do  i = 1, n
      if (dict%key(i) == key) then
         value = dict%value(i)
         status = 0
         RETURN
      endif
enddo

end subroutine wxml_get_value

!------------------------------------------------------
subroutine get_key(dict,i,key,status)
!
! Get the i'th key
!
type(wxml_dictionary_t), intent(in)            :: dict
integer, intent(in)                       :: i
character(len=*), intent(out)             :: key
integer, intent(out)                      :: status

if (i <= dict%number_of_items) then
      key = dict%key(i)
      status = 0
else
      key = ""
      status = -1
endif

end subroutine get_key

!------------------------------------------------------
subroutine add_key_to_dict(key,dict)
character(len=*), intent(in)          :: key
type(wxml_dictionary_t), intent(inout)   :: dict

integer  :: n

n = dict%number_of_items
if (n == MAX_ITEMS) then
      write(unit=0,fmt=*) "Dictionary capacity exceeded !"
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
character(len=*), intent(in)          :: value
type(wxml_dictionary_t), intent(inout)   :: dict

integer  :: n

n = dict%number_of_items
dict%value(n) = value

end subroutine add_value_to_dict

!------------------------------------------------------
subroutine reset_dict(dict)
type(wxml_dictionary_t), intent(inout)   :: dict

dict%number_of_items = 0

end subroutine reset_dict

!------------------------------------------------------
subroutine print_dict(dict)
type(wxml_dictionary_t), intent(in)   :: dict

integer  :: i

do i = 1, dict%number_of_items
      print *, trim(dict%key(i)), " = ", trim(dict%value(i))
enddo

end subroutine print_dict

end module xmlf90_wxml_dictionary
