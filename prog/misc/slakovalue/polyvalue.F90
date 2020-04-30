!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Reads a spline repulsive from an SK-table and returns its value and its first
!! and second derivatives.
program polyvalue
  use dftbp_accuracy
  use dftbp_io
  use dftbp_reppoly
  use dftbp_fileid
  use dftbp_message
  implicit none

  character(lc) :: arg, fname
  logical :: homo
  type(TRepPolyIn) :: repPolyIn
  type(TRepPoly) :: pRepPoly
  integer :: fp, iostat, ii, npoint
  real(dp), parameter :: rstart = 0.01_dp, dr = 0.01_dp
  real(dp) :: rr(3), energy, grad(3), d2, rDummy

  if (command_argument_count() == 0) then
    call error("Wrong number of arguments. Use 'polyvalue -h' to obtain help.")
  end if
  call get_command_argument(1, arg)
  if (arg == "-h" .or. arg == "--help") then
    write(stdout, "(A)") &
        & "Usage: polyvalue  homo | hetero  skfile",&
        & "",&
        & "Reads an SK-file, extracts the polynomial repulsive from it and &
        &prints its value", &
        & "and the first and second derivatives up to the repulsive cutoff. &
        &Output values",&
        & "are given in atomic units with Hartree as energy unit."
    stop
  end if
  if (arg /= "homo" .and. arg /= "hetero") then
    call error("The first argument must be 'homo' or 'hetero'")
  end if
  homo = (arg == "homo")
  if (command_argument_count() /= 2) then
    call error("Missing file name. Use 'polyvalue -h' to obtain help.")
  end if
  call get_command_argument(2, fname)

  fp = getFileId()
  open(fp, file=fname, action="read", status="old", iostat=iostat)
  if (iostat /= 0) then
    call error("Unable to open file '" // trim(fname) // "'")
  end if

  read(fp, *)
  if (homo) then
    read(fp, *)
  end if
  read(fp, *) rDummy, repPolyIn%polyCoeffs, repPolyIn%cutoff, &
      & (rDummy, ii = 11, 20)
  close(fp)

  call init(pRepPoly, repPolyIn)
  npoint = floor((repPolyIn%cutoff - rstart) / dr) + 1
  rr(:) = 0.0_dp
  do ii = 0, nPoint
    rr(1) = rStart + real(ii, dp) * dr
    call getenergy(pRepPoly, energy, rr(1))
    call getenergyderiv(pRepPoly, grad, rr, d2)
    write(stdout, "(4E23.15)") rr(1), energy, grad(1), d2
  end do

end program polyvalue
