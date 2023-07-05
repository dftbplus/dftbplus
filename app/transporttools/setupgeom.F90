!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2023  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program setupgeom
  use dftbp_common_globalenv
  use dftbp_io_formatout, only : printDftbHeader
  use transporttools_inputdata, only : TInputData
  use transporttools_parser, only : parseHsdInput
#:if WITH_MPI
  use mpi, only : MPI_THREAD_FUNNELED, MPI_COMM_WORLD
  use dftbp_common_mpienv, only : TMpiEnv, TMpiEnv_init
  use dftbp_extlibs_mpifx, only : mpifx_init_thread, mpifx_finalize
#:endif
  implicit none

  character(len=*), parameter :: releaseName = ''
  integer, parameter :: releaseYear = 2020

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
  call printDftbHeader(releaseName, releaseYear)
  allocate(input)
  call parseHsdInput(input)
  deallocate(input)
  call destructGlobalEnv()

#:if WITH_MPI
  call mpifx_finalize()
#:endif

end program setupgeom
