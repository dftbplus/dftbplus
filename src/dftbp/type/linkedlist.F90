!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains types and functions and subroutines for manipulating linked lists.  Every list must be
!> initialized with init, and destroyed with destroy.
module dftbp_type_linkedlist
  use dftbp_common_accuracy, only : lc, mc
  use dftbp_type_linkedlistc0, only : append, asArray, destruct, init, len, TListComplex
  use dftbp_type_linkedlistc1, only : append, asArray, asVector, destruct, init, intoArray, len,&
      & TListComplexR1
  use dftbp_type_linkedlisti0, only : append, asArray, destruct, init, len, TListInt
  use dftbp_type_linkedlisti1, only : append, asArray, asVector, destruct, elemShape, get, init,&
      & intoArray, len, TListIntR1
  use dftbp_type_linkedlistlc0, only : append, destruct, get, init, TListCharLc
  use dftbp_type_linkedlistmc0, only : TListCharMc
  use dftbp_type_linkedlistr0, only : append, asArray, destruct, init, len, TListReal
  use dftbp_type_linkedlistr1, only : append, asArray, asVector, destruct, init, intoArray, len,&
      & TListRealR1
  use dftbp_type_linkedlistr2, only : append, destruct, init, intoArray, len, TListRealR2
  use dftbp_type_linkedlists0, only : append, asArray, destruct, find, get, hasElement, init,&
      & isUnishaped, len, set, TListString
  implicit none

  private
  !> Expose the used linked list content
  public :: TListComplex, TListComplexR1
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
