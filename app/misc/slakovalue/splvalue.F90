!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2022  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Reads a spline repulsive from an SK-table and returns its value and its first
!! and second derivatives.
program splvalue
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_globalenv, only : stdOut
  use dftbp_dftb_repulsive_splinerep, only : TSplineRepInp, TSplineRep, TSplineRep_init
  use dftbp_io_fileid, only : getFileId
  use dftbp_io_message, only : error
  use dftbp_type_oldskdata, only : readsplinerep
#:if WITH_MPI
  use dftbp_common_mpienv, only : TMpiEnv
#:endif
  implicit none

  character(*), parameter :: fname = "test.skf"
  character(lc) :: arg
  type(TSplineRepInp) :: splineRepInp
  type(TSplineRep) :: splineRep
  integer :: fp, iostat, ii, npoint
  real(dp), parameter :: rstart = 0.01_dp, dr = 0.01_dp
  real(dp) :: rr, energy, dEnergy, d2Energy

#:if WITH_MPI
  !> MPI environment, if compiled with mpifort
  type(TMpiEnv) :: mpi

  ! As this is serial code, trap for run time execution on more than 1 processor with an mpi enabled
  ! build
  call TMpiEnv_init(mpi)
  call mpi%mpiSerialEnv()
#:endif

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
  call readsplinerep(fp, fname, splineRepInp)
  close(fp)

  call TSplineRep_init(splineRep, splineRepInp)
  npoint = floor((splineRepInp%cutoff - rstart) / dr) + 1
  do ii = 0, npoint
    rr = rstart + real(ii, dp) * dr
    call splineRep%getValue(rr, energy=energy, dEnergy=dEnergy, d2Energy=d2Energy)
    write(stdout, "(4E23.15)") rr, energy, dEnergy, d2Energy
  end do

end program splvalue
