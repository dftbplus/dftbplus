!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Inter-operation with C routines
module dftbp_common_clang
  use, intrinsic :: iso_c_binding, only : c_char, c_ptr, c_null_char, c_f_pointer, c_associated,&
      & c_size_t
  use, intrinsic :: iso_fortran_env, only : output_unit
  implicit none

  private
  public :: fortranChar, c_free, c_strlen

  !> Various C routines which need a Fortran binding
  interface

    !> Clean up memory attached to a C pointer
    subroutine c_free(ptr) bind(c, name="free")
      import :: c_ptr
      implicit none
      !> Pointer to nullify
      type(c_ptr), value :: ptr
    end subroutine c_free

    !> strlen for a dynamically allocated string
    integer(c_size_t) function c_strlen(s) bind(C, name='strlen')
      import :: c_size_t, c_ptr
      !> C string
      type(c_ptr), value :: s
    end function c_strlen

  end interface

contains

  !> Converts a 0-char terminated C-type string into a Fortran string.
  function fortranChar(cstring, maxlen)

    !> C-type string as array
    character(kind=c_char), intent(in) :: cstring(*)

    !> Maximal string length. If C-string is longer, it will be chopped.
    integer, intent(in), optional  :: maxlen

    !> Resulting Fortran string
    character(:, kind=c_char), allocatable :: fortranChar

    integer :: ii, maxlen0

    if (present(maxlen)) then
      maxlen0 = maxlen
    else
      maxlen0 = huge(maxlen0) - 1
    end if

    do ii = 1, maxlen0
      if (cstring(ii) == c_null_char) then
        exit
      end if
    end do
    allocate(character(ii - 1) :: fortranChar)
    fortranChar = transfer(cstring(1 : ii - 1), fortranChar)

  end function fortranChar


end module dftbp_common_clang
