!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Exporting mpifx functionality if compiled with mpi support, otherwise empty.
module dftbp_extlibs_mpifx
#:if WITH_MPI
  use libmpifx_module, only : MPI_COMM_TYPE_SHARED, MPI_LAND, MPI_LOR, MPI_MAX, MPI_SUM,&
      & MPI_THREAD_FUNNELED, mpifx_abort, mpifx_allgather, mpifx_allgatherv, mpifx_allreduce,&
      & mpifx_allreduceip, mpifx_barrier, mpifx_bcast, mpifx_comm, mpifx_finalize, mpifx_gatherv,&
      & mpifx_init_thread, mpifx_recv, mpifx_reduceip, mpifx_send, MPIFX_SIZE_T, mpifx_win
#:endif
  implicit none
  public

end module dftbp_extlibs_mpifx
