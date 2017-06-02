!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!!* Auxiliary subroutines for the ASSERT command
module assert
  implicit none
  private

  public :: assertError

contains  

  !!* Prints assertion error and abort program execution.
  !!* @param fileName Name of the file in which the error occured.
  !!* @param lineNr   Nr. of the line at which the error occured.
  !!* @desc
  !!*   Here follows the long description of the whole storry. Things like
  !!*   dependencies hasnt to be mentioned, because they should be generated
  !!*   automatically by the documentation system. Also things like type,
  !!*   intent etc. of the arguments should be generated automagically.
  !!*   (Nice dreams...).
  !!* @todo Maybe, the common error printing facility should be used.
  !!*       The documentation should be changed to RoboDoc style
  subroutine assertError(fileName, lineNr)
    character(*), intent(in) :: fileName
    integer,      intent(in) :: lineNr
    
    write (*, '(A)') "!!! UNFULLFILLED ASSERTION"
    write (*, '(A,A)') "!!! FILE:      ", fileName
    write (*, '(A,I0)') "!!! LINE NR.:  ", lineNr
    stop
    
  end subroutine assertError


end module assert
