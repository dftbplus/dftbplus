!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Module to print data types
module dftbp_io_intrinsicpr
  use dftbp_common_accuracy, only : dp, lc

  private
  public :: printContent


  !> Print various types of array
  interface printContent
    module procedure printArrayRealR1
    module procedure printArrayRealR2
    module procedure printArrayRealR3
    module procedure printArrayRealR4
    module procedure printArrayComplexR1
    module procedure printArrayComplexR2
    module procedure printArrayComplexR3
    module procedure printArrayComplexR4
    module procedure printArrayIntR1
    module procedure printArrayIntR2
    module procedure printArrayIntR3
    module procedure printArrayIntR4
    module procedure printArrayCharR1
    module procedure printArrayCharR2
  end interface printContent

contains

  ! Real arrays


  !> print real values
  subroutine printArrayRealR1(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    real(dp), intent(in) :: array(:)

    !> leave out header information
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    write(unit, *) (array(ii), ii = lbound(array, 1), ubound(array, 1))

  end subroutine printArrayRealR1


  !> print real values
  subroutine printArrayRealR2(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    real(dp), intent(in) :: array(:, :)

    !> leave out header information
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 2), ubound(array, 2)
      print *, "--2------", ii, "------"
      call printContent(unit, array(:, ii), .true.)
    end do

  end subroutine printArrayRealR2


  !> print real values
  subroutine printArrayRealR3(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    real(dp), intent(in) :: array(:, :, :)

    !> leave out header information
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 3), ubound(array, 3)
      print *, "--3------", ii, "------"
      call printContent(unit, array(:, :, ii), .true.)
    end do

  end subroutine printArrayRealR3


  !> print real values
  subroutine printArrayRealR4(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    real(dp), intent(in) :: array(:, :, :, :)

    !> leave out header information
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 4), ubound(array, 4)
      print *, "--4------", ii, "------"
      call printContent(unit, array(:, :, :, ii), .true.)
    end do

  end subroutine printArrayRealR4

  ! Complex arrays


  !> print complex values
  subroutine printArrayComplexR1(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    !> data to print
    complex(dp), intent(in) :: array(:)

    !> leave out header information
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    write(unit, *) (array(ii), ii = lbound(array, 1), ubound(array, 1))

  end subroutine printArrayComplexR1


  !> print complex values
  subroutine printArrayComplexR2(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    complex(dp), intent(in) :: array(:, :)

    !> leave out header information
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 2), ubound(array, 2)
      print *, "--2------", ii, "------"
      call printContent(unit, array(:, ii), .true.)
    end do

  end subroutine printArrayComplexR2


  !> print complex values
  subroutine printArrayComplexR3(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    complex(dp), intent(in) :: array(:, :, :)

    !> leave out header information
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 3), ubound(array, 3)
      print *, "--3------", ii, "------"
      call printContent(unit, array(:, :, ii), .true.)
    end do

  end subroutine printArrayComplexR3


  !> print complex values
  subroutine printArrayComplexR4(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    complex(dp), intent(in) :: array(:, :, :, :)

    !> leave out header information
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 4), ubound(array, 4)
      print *, "--4------", ii, "------"
      call printContent(unit, array(:, :, :, ii), .true.)
    end do

  end subroutine printArrayComplexR4

  ! Integer arrays


  !> print integer values
  subroutine printArrayIntR1(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    integer, intent(in) :: array(:)

    !> leave out header information
    logical, intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    write(unit, *) (array(ii), ii = lbound(array, 1), ubound(array, 1))

  end subroutine printArrayIntR1


  !> print integer values
  subroutine printArrayIntR2(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    integer, intent(in) :: array(:, :)

    !> leave out header information
    logical, intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 2), ubound(array, 2)
      print *, "--2------", ii, "------"
      call printContent(unit, array(:, ii), .true.)
    end do

  end subroutine printArrayIntR2


  !> print integer values
  subroutine printArrayIntR3(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    integer, intent(in) :: array(:, :, :)

    !> leave out header information
    logical, intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 3), ubound(array, 3)
      print *, "--3------", ii, "------"
      call printContent(unit, array(:, :, ii), .true.)
    end do

  end subroutine printArrayIntR3


  !> print integer values
  subroutine printArrayIntR4(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    integer, intent(in) :: array(:, :, :, :)

    !> leave out header information
    logical, intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 4), ubound(array, 4)
      print *, "--4------", ii, "------"
      call printContent(unit, array(:, :, :, ii), .true.)
    end do

  end subroutine printArrayIntR4

  ! Character arrays


  !> print character values
  subroutine printArrayCharR1(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    character(lc), intent(in) :: array(:)

    !> leave out header information
    logical,       intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    write(unit, *) (trim(array(ii)), ii = lbound(array, 1), ubound(array, 1))

  end subroutine printArrayCharR1


  !> print character values
  subroutine printArrayCharR2(unit, array, omitHeader)

    !> output unit
    integer, intent(in) :: unit

    !> data to print
    character(lc), intent(in) :: array(:, :)

    !> leave out header information
    logical,       intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 2), ubound(array, 2)
      print *, "--2------", ii, "------"
      call printContent(unit, array(:, ii), .true.)
    end do

  end subroutine printArrayCharR2

end module dftbp_io_intrinsicpr
