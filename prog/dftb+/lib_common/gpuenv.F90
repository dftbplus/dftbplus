!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Information on any GPUs on the system
module dftbp_gpuenv
  use iso_c_binding, only :  c_int
  use dftbp_magma, only : getGpusAvailable, getGpusRequested
  use dftbp_globalenv, only : stdOut
  implicit none

  private
  public :: TGpuEnv, TGpuEnv_init


  !> Contains global data for the GPU environment
  type :: TGpuEnv

    !> Number of GPUs available for the process
    integer(c_int) :: nGpu

  end type TGpuEnv


contains


  !> Initilises a TGpuEnv instance
  subroutine TGpuEnv_init(this)

    !> Instance
    type(TGpuEnv), intent(out) :: this

    integer(c_int) :: nGpuReq

    call getGpusAvailable(this%nGpu)
    call getGpusRequested(nGpuReq)
    write(stdOut, *) "Number of GPUs requested:", nGpuReq
    write(stdOut, *) "Number of GPUs found    :", this%nGpu
    if ((nGpuReq <= this%nGpu) .and. (nGpuReq >= 1)) then
      this%nGpu = nGpuReq
    end if

  end subroutine TGpuEnv_init


end module dftbp_gpuenv
