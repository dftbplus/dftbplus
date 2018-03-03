!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Module for filling in optional arguments if suplied to a call, otherwise suplying a default
!> value.
module optarg
  use accuracy, only : dp
  implicit none
  private

  public :: getOptionalArg


  !> Optional argument processor
  interface getOptionalArg
    module procedure getOptionalArgInt
    module procedure getOptionalArgReal
    module procedure getOptionalArgString
    module procedure getOptionalArgLogical
  end interface getOptionalArg

contains


  !> Return optional argument or default value if not present (int).
  subroutine getOptionalArgInt(defArg, outArg, optArg)

    !> Default value for optional argument.
    integer, intent(in) :: defArg

    !> Argument value on exit
    integer, intent(out) :: outArg

    !> Optional argument to check
    integer, intent(in), optional :: optArg

    if (present(optArg)) then
      outArg = optArg
    else
      outArg = defArg
    end if

  end subroutine getOptionalArgInt


  !> Return optional argument or default value if not present (real).
  subroutine getOptionalArgReal(defArg, outArg, optArg)

    !> Default value for optional argument.
    real(dp), intent(in) :: defArg

    !> Argument value on exit
    real(dp), intent(out) :: outArg

    !> Optional argument to check
    real(dp), intent(in), optional :: optArg

    if (present(optArg)) then
      outArg = optArg
    else
      outArg = defArg
    end if

  end subroutine getOptionalArgReal


  !> Return optional argument or default value if not present (str).
  subroutine getOptionalArgString(defArg, outArg, optArg)

    !> Default value for optional argument.
    character(*), intent(in) :: defArg

    !> Argument value on exit
    character(:), allocatable, intent(out) :: outArg

    !> Optional argument to check
    character(*), intent(in), optional :: optArg

    if (present(optArg)) then
      outArg = optArg
    else
      outArg = defArg
    end if

  end subroutine getOptionalArgString


  !> Return optional argument or default value if not present (logical).
  subroutine getOptionalArgLogical(defArg, outArg, optArg)

    !> Default value for optional argument.
    logical, intent(in) :: defArg

    !> Argument value on exit
    logical, intent(out) :: outArg

    !> Optional argument to check
    logical, intent(in), optional :: optArg

    if (present(optArg)) then
      outArg = optArg
    else
      outArg = defArg
    end if

  end subroutine getOptionalArgLogical

end module optarg
