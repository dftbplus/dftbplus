!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains types and functions and subroutines for manipulating linked lists.  Every list must be
!> initialized with init, and destroyed with destroy.
module linkedList
  use assert
  use accuracy
  use linkedlisti0
  use linkedlisti1
  use linkedlistr0
  use linkedlistr1
  use linkedlistr2
  use linkedlistmc0
  use linkedlistlc0
  use linkedlists0
  implicit none
  private


  !> Expose the used linked list content
  public :: listReal, listRealR1, listRealR2, listCharMc, listCharLc, listInt, listIntR1
  public :: listString
  public :: init, destruct, append, len, find, hasElement, elemShape, isUnishaped
  public :: get, set, asArray, asVector, intoArray
  public :: charMc, charLc

contains


  !> string of up to mc characters extraction from character array
  function charMc(string)

    !> full string
    character(*), intent(in) :: string

    !> resulting characters
    character(mc) :: charMC

    charMc = string(1:min(mc, len(string)))

  end function charMc


  !> string of up to lc characters extraction from character array
  function charLc(string)

    !> full string
    character(*), intent(in) :: string

    !> resulting characters
    character(lc) :: charLc

    charLc = string(1:min(lc, len(string)))

  end function charLc

end module linkedList
