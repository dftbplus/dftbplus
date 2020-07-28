!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Reads a spline repulsive from an SK-table and returns its value and its first
!! and second derivatives.
program integvalue
  use dftbp_accuracy
  use dftbp_globalenv, only : stdOut
  use dftbp_oldskdata
  use dftbp_slakoeqgrid
  use dftbp_fileid
  use dftbp_message
  implicit none

  integer, parameter :: nSKInter = 20
  integer, parameter :: nSKInterOld = 10
  integer, parameter :: iSKInterOld(nSKInterOld) &
      & = (/ 8, 9, 10, 13, 14, 15, 16, 18, 19, 20 /)
  real(dp), parameter :: deltaXDiff = epsilon(1.0_dp)**0.25_dp

  type(TOldSKData) :: skdata
  type(TSlakoEqGrid) :: skgrid
  character(lc) :: fname
  logical :: homo, extended
  integer :: nPoint, col
  real(dp), parameter :: rStart = 0.01_dp, dr = 0.001_dp
  real(dp), pointer :: data(:,:)

  call processArguments(fname, homo, extended, col)
  call readFromFile(skdata, fname, homo)
  call getSkColumnData(extended, skData, col, data)
  call init(skgrid, skdata%dist, data, skEqGridNew)
  nPoint = floor((getCutoff(skgrid) - rStart) / dr) + 1
  call writeValues(skgrid, rStart, dr, nPoint)

contains


  !> Prints help and stops.
  subroutine printHelp()

    write(stdout, "(A)") &
        & "Usage: integvalue  {homo|hetero}  skfile {orig|ext} col",&
        & "",&
        "Reads an SK-file, extracts the given column in the integral table and&
        & the ",&
        & "first and second derivatives up to the cutoff. Output values are&
        & given in ", &
        & "atomic units with Hartree as energy unit.", &
        & "",&
        & "homo|hetero -- whether SK file is homo or heteronuclear", &
        & "skfile      -- name of the SK file", &
        & "orig|ext    -- whether the column number is a meant in old or&
        & extended format", &
        & "col         -- column number"
    stop

  end subroutine printHelp


  !> Process program arguments.
  !!
  subroutine processArguments(fname, homo, extended, col)

    !> File name
    character(*), intent(out) :: fname

    !> homonuclear?
    logical, intent(out) :: homo

    !> extended format
    logical, intent(out) :: extended

    !> column to extract
    integer, intent(out) :: col

    character(lc) :: arg
    integer :: iostat

    if (command_argument_count() == 0) then
      call error("Wrong number of arguments. Use 'integvalue -h' to obtain&
          & help.")
    end if
    call get_command_argument(1, arg)
    if (arg == "-h" .or. arg == "--help") then
      call printHelp()
    end if
    if (command_argument_count() /= 4) then
      call error("Invalid number of arguments. Use 'integvalue -h' to obtain&
          & help.")
    end if
    if (arg /= "homo" .and. arg /= "hetero") then
      call error("The first argument must be 'homo' or 'hetero'")
    end if
    homo = (arg == "homo")

    call get_command_argument(2, fname)
    call get_command_argument(3, arg)
    if (arg /= "orig" .and. arg /= "ext" ) then
      call error("The third argument must be 'orig' or 'ext'")
    end if
    extended = (arg == "ext")
    call get_command_argument(4, arg)
    read(arg, *, iostat=iostat) col
    if (iostat /= 0) then
      call error("Third argument must the column number")
    end if

  end subroutine processArguments


  !> Returns the appropriate column of the SK-table.
  !!
  subroutine getSkColumnData(extended, skData, col, data)

    !> extended format
    logical, intent(in) :: extended

    !> Slater-Koster data
    type(TOldSKData), intent(in), target :: skData

    !> Column to extract
    integer, intent(in) :: col

    !> resulting data
    real(dp), pointer, intent(out) :: data(:,:)

    integer :: mycol

    if (extended) then
      if (col < 1 .or. col > 2 * nSKInter) then
        call error("Invalid column number")
      end if
      if (col <= nSkInter) then
        mycol = col
        data => skdata%skHam(:,mycol:mycol)
      else
        mycol = col - nSKInter
        data => skdata%skOver(:,mycol:mycol)
      end if
    else
      if (col < 1 .or. col > 2 * nSKInterOld) then
        call error("Invalid column number")
      end if
      if (col <= nSkInterOld) then
        mycol = iSKInterOld(col)
        data => skdata%skHam(:,mycol:mycol)
      else
        mycol = iSKInterOld(col - nSKInterOld)
        data => skdata%skOver(:,mycol:mycol)
      end if
    end if

  end subroutine getSkColumnData


  !> Writes values on a given grid.
  !!
  subroutine writeValues(skgrid, rStart, dr, nPoint)

    !> SK data grid
    type(TSlakoEqGrid), intent(in) :: skgrid

    !> starting distance
    real(dp), intent(in) :: rStart

    !> separation
    real(dp), intent(in) :: dr

    !> Number of points
    integer, intent(in) :: nPoint

    integer :: ii
    real(dp) :: dist, sk0(1), skm1(1), skp1(1)

    do ii = 0, nPoint
      dist = rStart + real(ii, dp) * dr
      call getSKIntegrals(skgrid, sk0, dist)
      call getSKIntegrals(skgrid, skp1, dist + deltaXDiff)
      call getSKIntegrals(skgrid, skm1, dist - deltaXDiff)
      write(stdout, "(4E23.15)") dist, sk0, (skp1 - skm1) / deltaXDiff, &
          & (skp1 + skm1 - 2.0_dp * sk0) / (deltaXDiff * deltaXDiff)
    end do

  end subroutine writeValues

end program integvalue
