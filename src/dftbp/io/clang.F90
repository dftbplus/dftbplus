!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Inter-operation with C routines
module dftbp_io_clang
  use, intrinsic :: iso_c_binding, only : c_char, c_ptr, c_null_char, c_f_pointer, c_associated
  use, intrinsic :: iso_fortran_env, only : output_unit
  implicit none

  private
  public :: fortranChar, handleOutputFileName

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


  !> Handles the optional output file name (which should be a NULL-ptr if not present), or use
  !> stdout instead
  subroutine handleOutputFileName(outputFileName, outputUnit, tOutputOpened)

    !> Name of file from C (or NULL)
    type(c_ptr), intent(in) :: outputFileName

    !> Resulting fortran unit
    integer, intent(out) :: outputUnit

    !> Was the file opened
    logical, intent(out) :: tOutputOpened

    character(c_char), pointer :: pOutputFileName
    character(:), allocatable :: fortranFileName

    if (c_associated(outputFileName)) then
      call c_f_pointer(outputFileName, pOutputFileName)
      fortranFileName = fortranChar(pOutputFileName)
      open(newunit=outputUnit, file=fortranFileName, action="write")
      tOutputOpened = .true.
    else
      outputUnit = output_unit
      tOutputOpened = .false.
    end if

  end subroutine handleOutputFileName

end module dftbp_io_clang
