!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Trivial MPI program to test compilable environment suports MPI
program test_mpiplatform
  use mpi
  implicit none

  integer, parameter :: requiredThreading = MPI_THREAD_FUNNELED
  integer :: providedThreading, rank, np, ierr
  integer, parameter :: lead_id = 0

  call mpi_init_thread(requiredThreading, providedThreading, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, np, ierr)

end program test_mpiplatform
