module xmlf90_charset
!
! One-byte only, sorry
!

  implicit none

private

integer, parameter, private  :: small_int = selected_int_kind(1)

!--------------------------------------------------------------------------
type, public :: charset_t
! private
      integer(kind=small_int), dimension(0:255) :: mask
end type charset_t


public  :: operator(.in.), operator(+)
public  :: assignment(=)
public  :: print_charset, reset_charset

interface operator(.in.)
      module procedure belongs
end interface
private :: belongs

interface assignment(=)
      module procedure set_string_to_charset, set_codes_to_charset
end interface
private ::  set_string_to_charset, set_codes_to_charset

interface operator(+)
      module procedure add_string_to_charset, &
                       add_code_to_charset, add_codes_to_charset
end interface
private :: add_string_to_charset, add_code_to_charset, add_codes_to_charset

!--------------------------------------------------------------------------

character(len=*), parameter, private :: &
            lowercase = "abcdefghijklmnopqrstuvwxyz",  & 
            uppercase = "ABCDEFGHIJKLMNOPQRSTUVWXYZ", &
            digits    = "0123456789"

integer, parameter, public   :: SPACE           = 32
integer, parameter, public   :: NEWLINE         = 10
integer, parameter, public   :: CARRIAGE_RETURN = 13
integer, parameter, public   :: TAB             =  9

type(charset_t), public    :: initial_name_chars
type(charset_t), public    :: name_chars
type(charset_t), public    :: whitespace
type(charset_t), public    :: valid_chars
type(charset_t), public    :: uppercase_chars

public  :: setup_xml_charsets


CONTAINS !==========================================================

!--------------------------------------------------------------
function belongs(c,charset)  result(res)
character(len=1), intent(in)  :: c
type(charset_t), intent(in)   :: charset
logical                       :: res

integer  :: code

code = ichar(c)
res = (charset%mask(code) == 1)

end function belongs

!--------------------------------------------------------------

function add_string_to_charset(charset,str)  result (sum)
type(charset_t), intent(in)      :: charset
character(len=*), intent(in)     :: str
type(charset_t)                  :: sum

integer :: length, code, i

sum%mask = charset%mask

length = len_trim(str)
do i = 1, length
      code = ichar(str(i:i))
      sum%mask(code) = 1
enddo
end function add_string_to_charset

!--------------------------------------------------------------

function add_code_to_charset(charset,code) result(sum)
type(charset_t), intent(in)      :: charset
integer, intent(in)              :: code
type(charset_t)                  :: sum

if ((code > 255) .or. (code < 0)) return
sum%mask = charset%mask
sum%mask(code) = 1

end function add_code_to_charset

!--------------------------------------------------------------
function add_codes_to_charset(charset,codes) result(sum)
type(charset_t), intent(in)        :: charset
integer, dimension(:), intent(in)  :: codes
type(charset_t)                    :: sum

integer  :: i

sum%mask = charset%mask
do i = 1, size(codes)
      if ((codes(i) > 255) .or. (codes(i) < 0)) cycle
      sum%mask(codes(i)) = 1
enddo
end function add_codes_to_charset

!--------------------------------------------------------------

subroutine set_string_to_charset(charset,str)      
type(charset_t), intent(out)   :: charset
character(len=*), intent(in)   :: str


integer :: length, code, i

charset%mask = 0

length = len_trim(str)
do i = 1, length
      code = ichar(str(i:i))
      charset%mask(code) = 1
enddo

end subroutine set_string_to_charset

!--------------------------------------------------------------

subroutine set_codes_to_charset(charset,codes)      
type(charset_t), intent(out)   :: charset
integer, dimension(:), intent(in)  :: codes

integer :: i

charset%mask = 0

do i = 1, size(codes)
      charset%mask(codes(i)) = 1
enddo

end subroutine set_codes_to_charset


!--------------------------------------------------------------
subroutine print_charset(charset)      
type(charset_t), intent(in)   :: charset

integer :: i

do i = 0, 255
      if (charset%mask(i) == 1) print *, "Code: ", i
enddo
end subroutine print_charset

!--------------------------------------------------------------

subroutine reset_charset(charset)      
type(charset_t), intent(inout)   :: charset

integer :: i

do i = 0, 255
      charset%mask(i) = 0
enddo
end subroutine reset_charset

!--------------------------------------------------------------

!--------------------------------------------------------
subroutine setup_xml_charsets()

integer :: i

uppercase_chars = uppercase
initial_name_chars = (lowercase  // uppercase //  "_:" )
name_chars = initial_name_chars + ( digits // ".-")
whitespace = (/ SPACE, NEWLINE, TAB, CARRIAGE_RETURN /)

valid_chars = whitespace + (/ (i, i=33,255) /)

end subroutine setup_xml_charsets
!--------------------------------------------------------

end module xmlf90_charset









