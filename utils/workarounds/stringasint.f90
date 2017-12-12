! PGI compilers interprets strings starting with with 'T' or 'F' as valid integers at read.
!
! Known to fail:
!   PGI 17.10
!
! Workaround: Check first character of the string, whether it can be an integer.
!
module mymod
  implicit none

contains

  subroutine testBug(str, iError)
    character(*), intent(in) :: str
    integer, intent(out) :: iError

    integer :: ii

    read(str, *, iostat=iError) ii

  end subroutine testBug


  subroutine testWorkaround(str, iError)
    character(*), intent(in) :: str
    integer, intent(out) :: iError

    integer :: ii

    if (validIntegerStart(str(1:1))) then
      read(str, *, iostat=iError) ii
    else
      iError = -1
    end if

  end subroutine testWorkaround


  !> Cheks whether a given character represents a valid staring character for an integer.
  pure function validIntegerStart(char) result(tValid)

    !> Character to check
    character, intent(in) :: char

    !> Whether it can be a starting character of a string representing an integer.
    logical :: tValid

    integer :: ind

    ind = iachar(char)
    ! Either a digit 0-9 or - or +
    tValid = ((ind >= 48 .and. ind <= 57) .or. ind == 45 .or. ind == 43)

  end function validIntegerStart

end module mymod


program test
  use mymod
  implicit none

  character(10) :: str
  integer :: iError

  str = "123"
  call testBug(str, iError)
  print *, "You should see zero (works):", iError
  str = "Telephone"
  call testBug(str, iError)
  print *, "You should see NON-ZERO (bug):", iError
  call testWorkaround(str, iError)
  print *, "You should see NON-ZERO (workaround):", iError

end program test
