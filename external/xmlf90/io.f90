module xmlf90_io

implicit none
  
!
! Basic  I/O tools
!
integer, public, save            :: io_eor, io_eof

public  :: get_unit, setup_io
private :: find_eor_eof

CONTAINS

! ----------------------------------------------------------------------
subroutine setup_io()
  call find_eor_eof(io_eor, io_eof)
end subroutine setup_io

! ----------------------------------------------------------------------
subroutine get_unit(lun,iostat)

! Get an available Fortran unit number

integer, intent(out)  :: lun
integer, intent(out)  :: iostat

integer :: i
logical :: unit_used

do i = 10, 99
   lun = i
   inquire(unit=lun,opened=unit_used)
   if (.not. unit_used) then
      iostat = 0
      return
   endif
enddo
iostat = -1
lun = -1
end subroutine get_unit
! ----------------------------------------------------------------------

subroutine find_eor_eof(io_eor,io_eof)
!
! Determines the values of the iostat values for End of File and 
! End of Record (in non-advancing I/O)
!
integer, intent(out)           :: io_eor
integer, intent(out)           :: io_eof

integer           :: lun, iostat
character(len=1)  :: c

call get_unit(lun,iostat)

if (iostat /= 0) stop "Out of unit numbers"

open(unit=lun,status="scratch",form="formatted", &
     action="readwrite",position="rewind",iostat=iostat)
if (iostat /= 0)   stop "Cannot open test file"

write(unit=lun,fmt=*)  "a"
write(unit=lun,fmt=*)  "b"

rewind(unit=lun)

io_eor = 0
do
  read(unit=lun,fmt="(a1)",advance="NO",iostat=io_eor) c
  if (io_eor /= 0) exit
enddo

io_eof = 0
do
  read(unit=lun,fmt=*,iostat=io_eof)
  if (io_eof /= 0) exit
enddo

!!!!!!!!print *, "IO_EOR, IO_EOF: ", io_eor, io_eof

close(unit=lun,status="delete")

end subroutine find_eor_eof

! ----------------------------------------------------------------------
end module xmlf90_io









