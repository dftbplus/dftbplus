module std_io
  implicit none
  private

  integer, public :: stdOut
  integer, public :: stdIn
  public ::  set_stdout, set_stdin
  
 contains

  subroutine set_stdout(stdo)
    integer, intent(in) :: stdo
    stdOut = stdo
  end subroutine set_stdout    

  subroutine set_stdin(stdI)
    integer, intent(in) :: stdI
    stdIn = stdI
  end subroutine set_stdin    



end module std_io
