!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2017  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module intrinsicpr
  use accuracy

  private

  public :: printContent


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
  end interface


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Real arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine printArrayRealR1(array, omitHeader)
    real(dp), intent(in)           :: array(:)
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    write (*, *) (array(ii), ii = lbound(array, 1), ubound(array, 1))

  end subroutine printArrayRealR1



  subroutine printArrayRealR2(array, omitHeader)
    real(dp), intent(in)           :: array(:, :)
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 2), ubound(array, 2)
      print *, "--2------", ii, "------"
      call printContent(array(:, ii), .true.)
    end do

  end subroutine printArrayRealR2



  subroutine printArrayRealR3(array, omitHeader)
    real(dp), intent(in)           :: array(:, :, :)
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 3), ubound(array, 3)
      print *, "--3------", ii, "------"
      call printContent(array(:, :, ii), .true.)
    end do

  end subroutine printArrayRealR3



  subroutine printArrayRealR4(array, omitHeader)
    real(dp), intent(in)           :: array(:, :, :, :)
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 4), ubound(array, 4)
      print *, "--4------", ii, "------"
      call printContent(array(:, :, :, ii), .true.)
    end do

  end subroutine printArrayRealR4



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Complex arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine printArrayComplexR1(array, omitHeader)
    complex(dp), intent(in)        :: array(:)
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    write (*, *) (array(ii), ii = lbound(array, 1), ubound(array, 1))

  end subroutine printArrayComplexR1



  subroutine printArrayComplexR2(array, omitHeader)
    complex(dp), intent(in)        :: array(:, :)
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 2), ubound(array, 2)
      print *, "--2------", ii, "------"
      call printContent(array(:, ii), .true.)
    end do

  end subroutine printArrayComplexR2



  subroutine printArrayComplexR3(array, omitHeader)
    complex(dp), intent(in)        :: array(:, :, :)
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 3), ubound(array, 3)
      print *, "--3------", ii, "------"
      call printContent(array(:, :, ii), .true.)
    end do

  end subroutine printArrayComplexR3



  subroutine printArrayComplexR4(array, omitHeader)
    complex(dp), intent(in)        :: array(:, :, :, :)
    logical,  intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 4), ubound(array, 4)
      print *, "--4------", ii, "------"
      call printContent(array(:, :, :, ii), .true.)
    end do

  end subroutine printArrayComplexR4



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Integer arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine printArrayIntR1(array, omitHeader)
    integer, intent(in)           :: array(:)
    logical, intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    write (*, *) (array(ii), ii = lbound(array, 1), ubound(array, 1))

  end subroutine printArrayIntR1



  subroutine printArrayIntR2(array, omitHeader)
    integer, intent(in)           :: array(:, :)
    logical, intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 2), ubound(array, 2)
      print *, "--2------", ii, "------"
      call printContent(array(:, ii), .true.)
    end do

  end subroutine printArrayIntR2



  subroutine printArrayIntR3(array, omitHeader)
    integer, intent(in)           :: array(:, :, :)
    logical, intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 3), ubound(array, 3)
      print *, "--3------", ii, "------"
      call printContent(array(:, :, ii), .true.)
    end do

  end subroutine printArrayIntR3



  subroutine printArrayIntR4(array, omitHeader)
    integer, intent(in)           :: array(:, :, :, :)
    logical, intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 4), ubound(array, 4)
      print *, "--4------", ii, "------"
      call printContent(array(:, :, :, ii), .true.)
    end do

  end subroutine printArrayIntR4



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   Character arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine printArrayCharR1(array, omitHeader)
    character(lc), intent(in)           :: array(:)
    logical,       intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    write (*, *) (trim(array(ii)), ii = lbound(array, 1), ubound(array, 1))

  end subroutine printArrayCharR1



  subroutine printArrayCharR2(array, omitHeader)
    character(lc), intent(in)           :: array(:, :)
    logical,       intent(in), optional :: omitHeader

    integer :: ii

    if (.not. present(omitHeader)) then
      print *, " Shape: ", shape(array)
    end if
    do ii = lbound(array, 2), ubound(array, 2)
      print *, "--2------", ii, "------"
      call printContent(array(:, ii), .true.)
    end do

  end subroutine printArrayCharR2


end module intrinsicpr

