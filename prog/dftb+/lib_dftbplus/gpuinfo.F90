!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'


!> Information on any GPUs on the system
module dftbp_gpuinfo

#:if WITH_GPU
  use iso_c_binding, only :  c_int
  use device_info
  use dftbp_globalenv, only : stdOut
#:else
  use dftbp_message, only : error
#:endif
  implicit none

#:if WITH_GPU
  integer (c_int), protected :: ngpus
  integer (c_int), protected :: req_ngpus
#:endif

contains
  
  subroutine gpuInfo()

#:if WITH_GPU
    call  gpu_avail(ngpus)
    call  gpu_req(req_ngpus)
    write(stdOut,*) "Number of GPUs requested:",req_ngpus
    write(stdOut,*) "Number of GPUs found    :",ngpus
    if ((req_ngpus .le. ngpus) .and. (req_ngpus .ge. 1)) then
      ngpus = req_ngpus
    endif
#:else
    call error("Compiled without GPU support")
#:endif

  end subroutine gpuInfo


end module dftbp_gpuinfo
