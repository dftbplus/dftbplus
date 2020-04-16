!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Reads a spline repulsive from an SK-table and returns its value and its first
!! and second derivatives.
program splvalue
  use dftbp_accuracy
  use dftbp_io
  use dftbp_repspline
  use dftbp_oldskdata, only : readsplinerep
  use dftbp_fileid
  use dftbp_message
  implicit none

  character(*), parameter :: fname = "test.skf"
  character(lc) :: arg
  type(TRepSplinein) :: repsplinein
  type(TRepSpline) :: prepspline
  integer :: fp, iostat, ii, npoint
  real(dp), parameter :: rstart = 0.01_dp, dr = 0.01_dp
  real(dp) :: rr(3), energy, grad(3), d2

  if (command_argument_count() /= 1) then
    call error("Wrong number of arguments. Use 'splvalue -h' to obtain help.")
  end if
  call get_command_argument(1, arg)
  if (arg == "-h" .or. arg == "--help") then
    write(stdout, "(A)") &
        & "Usage: splvalue [ options ] skfile",&
        & "",&
        & "Reads an SK-file, extracts the spline repulsive from it and prints &
        &its value", &
        & "and the first and second derivatives up to the repulsive cutoff. &
        &Output values",&
        & "are given in atomic units with Hartree as energy unit."
    stop
  end if

  fp = getfileid()
  open(fp, file=arg, action="read", status="old", iostat=iostat)
  if (iostat /= 0) then
    call error("Unable to open file '" // trim(fname) // "'")
  end if
  call readsplinerep(fp, fname, repsplinein)
  close(fp)

  call init(prepspline, repsplinein)
  npoint = floor((repsplinein%cutoff - rstart) / dr) + 1
  rr(:) = 0.0_dp
  do ii = 0, npoint
    rr(1) = rstart + real(ii, dp) * dr
    call getenergy(prepspline, energy, rr(1))
    call getenergyderiv(prepspline, grad, rr, d2)
    write(stdout, "(4E23.15)") rr(1), energy, grad(1), d2
  end do

end program splvalue
