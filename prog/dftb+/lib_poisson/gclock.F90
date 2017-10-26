module gclock

  implicit none
  private
  
 integer, PARAMETER :: NMAXCLKS=5
 integer, public, save :: t1(NMAXCLKS),t2(NMAXCLKS),cr,cm,nclks=0

 public :: message_clock, set_clock, write_clock
 
contains

 subroutine message_clock(message)

   character(*) :: message
   integer :: l_mess
   character(2) :: str_mess, str_dots

   l_mess=len(message)
   write(str_mess,'(I2)') l_mess
   write(str_dots,'(I2)') 54-l_mess

   write(*,FMT='(A'//str_mess//','//str_dots//'("."))',ADVANCE='NO') message 
   
   !call flush(6)

   call set_clock

 end subroutine message_clock

 subroutine set_clock

   if (nclks.lt.NMAXCLKS) then
      nclks = nclks + 1
      call SYSTEM_CLOCK(t1(nclks),cr,cm) 
   endif

 end subroutine set_clock


 subroutine write_clock

   if (nclks.gt.0) then
      call SYSTEM_CLOCK(t2(nclks),cr,cm) 
      write(*,*) (t2(nclks)-t1(nclks))*1.0/cr,"sec"
      nclks=nclks-1
   endif

 end subroutine write_clock


 


end module gclock
