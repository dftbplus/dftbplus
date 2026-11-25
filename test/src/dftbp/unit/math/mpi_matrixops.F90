!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include "fortuno_mpi.fypp"

module test_math_mpi_matrixops
  use mpi_f08
  use fortuno_mpi, only : suite => mpi_suite_item, mpi_context, test_list, num_ranks, this_rank,&
      & global_comm
  $:FORTUNO_MPI_IMPORTS()
  implicit none (type, external)
  private

  public :: tests

contains

  !> Test an MPI_Allreduce sum across all ranks in the test communicator
  $:TEST("allreduce_sum")

    integer :: local_val, global_sum, expected_sum, nprocs, ierr

    ! Each rank contributes (rank + 1)
    local_val = this_rank() + 1
    nprocs = num_ranks()

    call MPI_Allreduce(local_val, global_sum, 1, MPI_INTEGER, MPI_SUM, global_comm(), ierr)

    ! Sum of 1..N is N*(N+1)/2
    expected_sum = nprocs * (nprocs + 1) / 2

    @:ASSERT(ierr == 0)
    @:ASSERT(global_sum == expected_sum)

  $:END_TEST()


  !> Test an MPI_Bcast from root rank 0
  $:TEST("bcast_value")

    integer :: val, ierr

    if (this_rank() == 0) then
      val = 42
    else
      val = -1
    end if

    call MPI_Bcast(val, 1, MPI_INTEGER, 0, global_comm(), ierr)

    @:ASSERT(ierr == 0)
    @:ASSERT(val == 42)

  $:END_TEST()


  !> Automatically collect all $:TEST(...) items in this module
  function tests()
    type(test_list) :: tests

    tests = test_list([&
        suite("mpi_matrixops", test_list([&
            $:TEST_ITEMS()
        ]))&
    ])
    $:STOP_ON_MISSING_TEST_ITEMS()

  end function tests

end module test_math_mpi_matrixops
