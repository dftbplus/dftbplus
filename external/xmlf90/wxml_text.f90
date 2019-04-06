module xmlf90_wxml_text

  implicit none

!
integer, private, parameter ::  sp = selected_real_kind(6,30)
integer, private, parameter ::  dp = selected_real_kind(14,100)
!
! TODO : Add optional format parameter
!

private

public :: str

interface str
   module procedure  str_integer,  str_real_dp, str_real_sp, &
                     str_logical
end interface
private :: str_integer,  str_real_dp, str_real_sp, str_logical

CONTAINS

      function str_integer(int,format) result(s)
      integer, intent(in)   :: int
      character(len=*), intent(in), optional  :: format
      character(len=100)    :: s

      if (present(format)) then
         write(s,format) int
      else
         write(s,"(i25)") int
      endif
      s = adjustl(s)
      end function str_integer

      function str_logical(log,format) result(s)
      logical, intent(in)   :: log
      character(len=*), intent(in), optional  :: format
      character(len=100)    :: s

      if (present(format)) then
         write(s,format) log
      else
         write(s,"(l1)") log
      endif
      s = adjustl(s)
      end function str_logical

      function str_real_dp(x,format) result(s)
      real(kind=dp), intent(in)   :: x
      character(len=*), intent(in), optional  :: format
      character(len=100)    :: s

      if (present(format)) then
         write(s,format) x
      else
         write(s,"(g22.12)") x
      endif
      s = adjustl(s)
      end function str_real_dp

      function str_real_sp(x,format) result(s)
      real(kind=sp), intent(in)   :: x
      character(len=*), intent(in), optional  :: format
      character(len=100)    :: s

      if (present(format)) then
         write(s,format) x
      else
         write(s,"(g22.12)") x
      endif
      s = adjustl(s)
      end function str_real_sp


end module xmlf90_wxml_text
