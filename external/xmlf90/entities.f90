module xmlf90_entities
  
!
! Entity management
!
! It deals with: 
!    1. The five standard entities (gt,lt,amp,apos,quot)
!    2. Character entities  (but only within the range of the char intrinsic)
!
use xmlf90_buffer

  implicit none

private

integer, parameter, private      :: MAX_REPLACEMENT_SIZE = 200
!
type, private :: entity_t
 character(len=40)                     :: code
 character(len=MAX_REPLACEMENT_SIZE)   :: replacement
end type entity_t

integer, parameter, private                          ::  N_ENTITIES  = 5

type(entity_t), private, dimension(N_ENTITIES), save :: predefined_ent =  &
      (/                  &
      entity_t("gt",">"), &
      entity_t("lt","<"),  &
      entity_t("amp","&"),  &
      entity_t("apos","'"),  &
      entity_t("quot","""")  &
      /)
      
public :: code_to_str , entity_filter

CONTAINS

subroutine code_to_str(code,str,status)
character(len=*), intent(in)  :: code
character(len=*), intent(out) :: str
integer, intent(out)          :: status         
integer :: i

integer   :: number, ll
character(len=4)  :: fmtstr

status = -1
do i = 1, N_ENTITIES
      if (code == predefined_ent(i)%code) then
         str = predefined_ent(i)%replacement
         status = 0
         return
      endif
enddo
!
! Replace character references  (but only within the range of the
! char intrinsic !!)
!
if (code(1:1) == "#") then
   if (code(2:2) == "x") then       ! hex character reference
      ll = len_trim(code(3:))
      write(unit=fmtstr,fmt="(a2,i1,a1)") "(Z", ll,")"
      read(unit=code(3:),fmt=fmtstr) number
      str = char(number)
      status = 0
      return
   else                             ! decimal character reference
      read(unit=code(2:),fmt=*) number
      str = char(number)
      status = 0
      return
   endif
endif

end subroutine code_to_str

!----------------------------------------------------------------
!
! Replaces entity references in buf1 and creates a new buffer buf2.
!
subroutine entity_filter(buf1,buf2,status,message)
type(buffer_t), intent(in)    :: buf1
type(buffer_t), intent(out)   :: buf2
integer, intent(out)          :: status
character(len=*), intent(out) :: message
!
! Replaces entity references by their value
!
integer :: i, k, len1
character(len=MAX_BUFF_SIZE)           :: s1
character(len=1)                       :: c
character(len=MAX_REPLACEMENT_SIZE)    :: repl

call buffer_to_character(buf1,s1)        !! Avoid allocation of temporary
len1 = len(buf1)

i = 1
status = 0

call reset_buffer(buf2)

do
   if (i > len1) exit
   c = s1(i:i)
   if (c == "&") then
      if (i+1 > len1) then
         status = -i
         message=  " Unmatched & in entity reference"
         return
      endif
      k = index(s1(i+1:),";")
      if (k == 0) then
         status = -i
         message=  " Unmatched & in entity reference"
         return
      endif
      call code_to_str(s1(i+1:i+k-1),repl,status)
      if (status /= 0) then
         status =  i     ! Could let it continue
         message= "Ignored unknown entity: &" // s1(i+1:i+k-1) // ";"
      else
         call add_to_buffer(trim(repl),buf2)
      endif
      i = i + k + 1
   else
     call add_to_buffer(c,buf2)
     i = i + 1
   endif
enddo

end subroutine entity_filter

end module xmlf90_entities





