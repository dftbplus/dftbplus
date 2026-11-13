!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

! Splits MPI_COMM_WORLD into two groups (round-robin, so that the ranks of the group members differ
! from their ranks in MPI_COMM_WORLD). Only the first group calculates H2O and Si2 systems via DFTB+
! (creating and destroying the DFTB+ instance in each cycle), while the processes of the second
! group wait on their own communicator meanwhile.
!
! Test whether DFTB+ restricts all its MPI and BLACS communication to the MPI communicator passed
! via the API. Any collective operation involving processes outside of that communicator (e.g. one
! on MPI_COMM_WORLD or on a BLACS grid derived from the default system context) would hang.

program test_mpisubgrids_idle
  ! No explicit import list, as some mpi modules miss symbols
  use mpi
  use, intrinsic :: iso_fortran_env, only : real64
  use dftbplus, only : TDftbPlus, TDftbPlus_init, TDftbPlusInput
  use testhelpers, only : writeAutotestTag   ! Only needed for the internal test system
  use dftbp_common_constants, only : AA__Bohr   ! Imported to ensure accurate conversion
  ! No explicit import list, as some mpi modules miss symbols
  implicit none

  integer, parameter :: dp = real64

  real(dp), parameter :: coordsSi2(3, 2) = reshape([&
      &  0.00000000000E+00_dp,  0.00000000000E+00_dp,  0.00000000000E+00_dp,&
      &  0.14567730000E+01_dp,  0.10567730000E+01_dp,  0.19567730000E+01_dp], [3, 2]) * AA__Bohr

  real(dp), parameter :: latVecsSi2(3, 3) = reshape([&
      &  0.22135460000E+01_dp,  0.24135460000E+01_dp,  0.00000000000E+00_dp,&
      &  0.00000000000E+00_dp,  0.27135460000E+01_dp,  0.27135460000E+01_dp,&
      &  0.29135460000E+01_dp,  0.00000000000E+00_dp,  0.20135460000E+01_dp], [3, 3]) * AA__Bohr

  real(dp), parameter :: coordsH2O(3, 3) = reshape([&
      & 0.00000000000E+00_dp,  -0.10000000000E+01_dp,   0.00000000000E+00_dp,&
      & 0.00000000000E+00_dp,   0.00000000000E+00_dp,   0.88306400000E+00_dp,&
      & 0.00000000000E+00_dp,   0.00000000000E+00_dp,  -0.78306400000E+00_dp], [3, 3]) * AA__Bohr

  !> Number of calculations (must be even, so that H2O and Si2 are calculated equally often)
  integer, parameter :: nRepeat = 10

  integer, parameter :: requiredThreading = MPI_THREAD_FUNNELED

  !> Number of groups to create
  integer, parameter :: nGroup = 2

  !> Group doing the DFTB+ calculations
  integer, parameter :: workingGroup = 0

  call main_()

contains


  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitely
  !! deallocated before the program finishes  (avoiding residual memory).
  !!
  subroutine main_()

    type(TDftbPlus), allocatable :: dftbp
    type(TDftbPlusInput), allocatable :: input

    real(dp) :: merminEnergy, totalMerminEnergy
    real(dp) :: gradients(3, 3), totalGradients(3, 3)
    integer :: devNull
    integer :: providedThreading, iErr, myId, myGroup, nProc
    integer :: iRepeat, nAtom
    integer :: groupComm, myGroupId, nGroupProc
    logical :: doSi2

    open(newunit=devNull, file="/dev/null", action="write")

    call mpi_init_thread(requiredThreading, providedThreading, iErr)
    call mpi_comm_rank(MPI_COMM_WORLD, myId, iErr)
    call mpi_comm_size(MPI_COMM_WORLD, nProc, iErr)

    if (modulo(nProc, nGroup) /= 0) then
      error stop "Number of processors must be a multiple of the number of groups"
    end if

    ! Processes are assigned to group in round-Robin fashion
    myGroup = modulo(myId, nGroup)
    call mpi_comm_split(MPI_COMM_WORLD, myGroup, myId / nGroup, groupComm, iErr)
    call mpi_comm_rank(groupComm, myGroupId, iErr)
    call mpi_comm_size(groupComm, nGroupProc, iErr)

    print "('ID: ', I2.2, ' | Group: ', I2.2, ' | ID in group: ', I2.2, ' (of ', I2.2, ')')",&
        & myId, myGroup, myGroupId, nGroupProc

    totalMerminEnergy = 0.0_dp
    totalGradients(:,:) = 0.0_dp

    if (myGroup == workingGroup) then

      repeat: do iRepeat = 1, nRepeat

        allocate(dftbp)
        allocate(input)

        doSi2 = modulo(iRepeat, 2) == 0
        if (doSi2) then
          nAtom = 2
        else
          nAtom = 3
        end if

        call TDftbPlus_init(dftbp, mpiComm=groupComm, devNull=devNull)

        if (doSi2) then
          call dftbp%getInputFromFile("dftb_in.Si2.hsd", input)
        else
          call dftbp%getInputFromFile("dftb_in.H2O.hsd", input)
        end if

        call dftbp%setupCalculator(input)

        if (doSi2) then
          call dftbp%setGeometry(coordsSi2, latVecsSi2)
        else
          call dftbp%setGeometry(coordsH2O)
        end if

        call dftbp%getEnergy(merminEnergy)
        totalMerminEnergy = totalMerminEnergy + merminEnergy
        print "('[', I2.2, '|', I2.2, '] (', I3.3, ') Obtained Mermin Energy:', F15.10)", myGroup,&
            & myGroupId, iRepeat, merminEnergy

        call dftbp%getGradients(gradients(:, :nAtom))
        totalGradients(:, :nAtom) = totalGradients(:, :nAtom) + gradients(:, :nAtom)

        ! Deallocating instances to trigger finalizers
        deallocate(dftbp)
        deallocate(input)

      end do repeat

    end if

    ! Processes of the idle group wait here (on their own communicator) while the working group
    ! does the DFTB+ calculations.
    call mpi_barrier(groupComm, iErr)

    ! The lead of the working group writes the averaged results
    if (myGroup == workingGroup .and. myGroupId == 0) then
      totalMerminEnergy = totalMerminEnergy / real(nRepeat, dp)
      totalGradients(:,:) = totalGradients / real(nRepeat, dp)
      print *, "AVERAGE MERMIN:"
      print *, totalMerminEnergy
      print *, "AVERAGE GRAD:"
      print *, totalGradients

      ! Write file for internal test system
      call writeAutotestTag(merminEnergy=totalMerminEnergy, gradients=totalGradients)
    end if

    call mpi_barrier(MPI_COMM_WORLD, iErr)
    call mpi_comm_free(groupComm, iErr)
    call mpi_finalize(iErr)

  end subroutine main_

end program test_mpisubgrids_idle
