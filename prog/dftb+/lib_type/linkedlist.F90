!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains types and functions and subroutines for manipulating linked lists.  Every list must be
!> initialized with init, and destroyed with destroy.
module dftbp_type_linkedlist
  use dftbp_common_assert
  use dftbp_common_accuracy
  use dftbp_type_linkedlisti0
  use dftbp_type_linkedlisti1
  use dftbp_type_linkedlistr0
  use dftbp_type_linkedlistr1
  use dftbp_type_linkedlistr2
  use dftbp_type_linkedlistmc0
  use dftbp_type_linkedlistlc0
  use dftbp_type_linkedlists0
  implicit none
  private


  !> Expose the used linked list content
  public :: TListReal, TListRealR1, TListRealR2
  public :: TListCharMc, TListCharLc, TListInt, TListIntR1
  public :: TListString
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

end module dftbp_type_linkedlist
