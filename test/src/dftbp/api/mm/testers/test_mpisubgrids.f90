!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

! Splits MPI_COMM_WORLD into groups of four processors. Each group calculates repeatedly H2O and
! Si2 systems via DFTB+, creating and destroying the DFTB+ instance in each cycle.
! This should test the capability of running multiple DFTB+ distances via the API on
! non-overlapping MPI subgrids.

program test_mpisubgrids
  ! Omitting explicit import list, as some mpi modules simply miss symbols (as in MPICH 3.3)
  use mpi
      !, only : MPI_COMM_WORLD, MPI_THREAD_FUNNELED, MPI_REAL8, MPI_SUM, mpi_comm,&
      !& mpi_init_thread, mpi_comm_rank, mpi_comm_size, mpi_comm_split, mpi_ireduce, mpi_barrier,&
      !& mpi_comm_free, mpi_finalize
  use dftbplus, only : TDftbPlusInput, TDftbPlus, TDftbPlus_init
  use dftbp_common_constants, only : AA__Bohr   ! Imported to ensure accurate conversion
  use testhelpers, only : writeAutotestTag   ! Only needed for the internal test system
  use, intrinsic :: iso_fortran_env, only : real64
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

  integer, parameter :: nRepeat = 100

  integer, parameter :: requiredThreading = MPI_THREAD_FUNNELED

  integer, parameter :: groupSize = 4

  call main_()

contains


  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitely
  !! deallocated before the program finishes.
  !!
  subroutine main_()

    type(TDftbPlus), allocatable :: dftbp
    type(TDftbPlusInput), allocatable :: input

    real(dp) :: merminEnergy, totalMerminEnergy, reducedMerminEnergy
    real(dp) :: gradients(3, 3), totalGradients(3, 3), reducedGradients(3, 3)
    integer :: devNull
    integer :: providedThreading, iErr, myId, myGroup, nProc
    integer :: iRepeat, nAtom

    !type(mpi_comm) :: groupComm
    integer :: groupComm
    integer :: nGroup, myIdGroup, nProcGroup

    logical :: doSi2

    open(newunit=devNull, file="/dev/null", action="write")

    call mpi_init_thread(requiredThreading, providedThreading, iErr)
    call mpi_comm_rank(MPI_COMM_WORLD, myId, iErr)
    call mpi_comm_size(MPI_COMM_WORLD, nProc, iErr)

    if (modulo(nProc, groupSize) /= 0) then
      error stop "Number of processors must be a multiple of the groupSize"
    end if

    nGroup = nProc / groupSize
    if (myId == 0) then
      print "('Creating ', I0, ' group(s)')", nGroup
    end if

    ! Processes are assigned to group in round-Robin fashion
    myGroup = modulo(myId, nGroup)
    call mpi_comm_split(MPI_COMM_WORLD, modulo(myId, nGroup), myId / nGroup, groupComm, iErr)
    call mpi_comm_rank(groupComm, myIdGroup, iErr)
    call mpi_comm_size(groupComm, nProcGroup, iErr)

    print "('ID: ', I2.2, ' | Group: ', I2.2, ' | ID in group: ', I2.2, ' (of ', I2.2, ')')",&
        & myId, myGroup, myIdGroup, nProcGroup

    totalMerminEnergy = 0.0_dp
    totalGradients(:,:) = 0.0_dp

    repeat: do iRepeat = 1, nRepeat

      allocate(dftbp)
      allocate(input)

      if (myId == 0) then
        print "('*** Cycle ', I0)", iRepeat
      end if

      doSi2 = modulo(iRepeat + myGroup, 2) == 0
      if (doSi2) then
        nAtom = 2
      else
        nAtom =3
      end if

      call TDftbPlus_init(dftbp, mpiComm=groupComm, devNull=devNull)
      !call TDftbPlus_init(dftbp, outputUnit=devNull, mpiComm=groupComm, devNull=devNull)

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
          & myIdGroup, iRepeat, merminEnergy

      call dftbp%getGradients(gradients(:, :nAtom))
      totalGradients(:, :nAtom) = totalGradients(:, :nAtom) + gradients(:, :nAtom)
      print "('[', I2.2, '|', I2.2, '] (', I3.3, ') Obtained gradient of atom 1:', 3F15.10)",&
          & myGroup, myIdGroup, iRepeat, gradients(:, 1)

      ! Deallocating instances to trigger finalizers
      deallocate(dftbp)
      deallocate(input)

    end do repeat

    call mpi_reduce(totalMerminEnergy, reducedMerminEnergy, 1, MPI_REAL8, MPI_SUM, 0,&
        & MPI_COMM_WORLD, iErr)
    call mpi_reduce(totalGradients, reducedGradients, size(totalGradients), MPI_REAL8, MPI_SUM, 0,&
        & MPI_COMM_WORLD, iErr)

    if (myId == 0) then
      reducedMerminEnergy = reducedMerminEnergy / real(nProc * nRepeat, dp)
      reducedGradients(:,:) = reducedGradients / real(nProc * nRepeat, dp)
      print *, "REDUCED MERMIN:"
      print *, reducedMerminEnergy
      print *, "REDUCED GRAD:"
      print *, reducedGradients

      ! Write file for internal test system
      call writeAutotestTag(merminEnergy=reducedMerminEnergy, gradients=reducedGradients)
    end if

    call mpi_barrier(MPI_COMM_WORLD, iErr)
    call mpi_comm_free(groupComm, iErr)
    call mpi_finalize(iErr)

  end subroutine main_

end program test_mpisubgrids

! H2O
! Total Mermin free energy:           -3.9711540041 H         -108.0606 eV
! Total Forces
!          0.000000000000      0.180238226181      0.008042274776
!          0.000000000000     -0.089537143031     -0.052333465531
!         -0.000000000000     -0.090701083150      0.044291190756


! Si2
! Total Mermin free energy:           -2.4599532761 H          -66.9387 eV
! Total Forces
!          0.089808060651     -0.050576358492      0.032427025640
!         -0.089808060651      0.050576358492     -0.032427025640
!

! Average:
! Mermin free: -3.2155536401
! Total GRADIENTS
!   -0.0449040303255 -0.0648309338445      -0.020234650208000002
!    0.0449040303255  0.019480392269499998  0.0423802455855
!    0.0              0.045350541575       -0.022145595378
