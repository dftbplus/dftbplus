!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains types and functions and subroutines for manipulating linked lists.  Every list must be
!> initialized with init, and destroyed with destroy.
module dftbp_type_linkedlist
  use dftbp_common_accuracy, only : mc, lc
  use dftbp_type_linkedlistc0, only : TListComplex, append, init, len, destruct, asArray
  use dftbp_type_linkedlistc1, only : TListComplexR1, append, len, init, asArray, destruct,&
      & asVector, intoArray
  use dftbp_type_linkedlisti0, only : TListInt, append, destruct, asArray, init, len
  use dftbp_type_linkedlisti1, only : TListIntR1, append, len, destruct, init, asArray, asVector,&
      & intoArray, get, elemShape
  use dftbp_type_linkedlistlc0, only : TListCharLc, append, init, get, destruct
  use dftbp_type_linkedlistmc0, only : TListCharMc
  use dftbp_type_linkedlistr0, only : TListReal, append, init, len, destruct, asArray
  use dftbp_type_linkedlistr1, only : TListRealR1, append, len, init, asArray, destruct, asVector,&
      & intoArray
  use dftbp_type_linkedlistr2, only : TListRealR2, init, append, intoArray, len, destruct
  use dftbp_type_linkedlists0, only : TListString, find, set, hasElement, isUnishaped, append, len,&
      & init, asArray, destruct, get
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
