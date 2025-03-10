!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program setupgeom
  use transporttools_inputdata, only : TInputData
  use transporttools_parser, only : parseHsdInput
  use dftbp_common_globalenv, only : destructGlobalEnv, initGlobalEnv
  use dftbp_common_release, only : releaseYear
  use dftbp_io_formatout, only : printDftbHeader
#:if WITH_MPI
  use mpi, only : MPI_COMM_WORLD, MPI_THREAD_FUNNELED
  use dftbp_common_mpienv, only : TMpiEnv, TMpiEnv_init
  use dftbp_extlibs_mpifx, only : mpifx_finalize, mpifx_init_thread
#:endif
  implicit none

  type(TInputData), allocatable :: input

#:if WITH_MPI
  !> MPI environment, if compiled with mpifort
  type(TMpiEnv) :: mpiEnv

  ! As this is serial code, trap for run time execution on more than 1 processor with an mpi enabled
  ! build
  call mpifx_init_thread(requiredThreading=MPI_THREAD_FUNNELED)
  call TMpiEnv_init(mpiEnv)
  call mpiEnv%mpiSerialEnv()
  call initGlobalEnv(mpiComm=MPI_COMM_WORLD)
#:else
  call initGlobalEnv()
#:endif
  call printDftbHeader('(setupgeom)', releaseYear)
  allocate(input)
  call parseHsdInput(input)
  deallocate(input)
  call destructGlobalEnv()

#:if WITH_MPI
  call mpifx_finalize()
#:endif

end program setupgeom
