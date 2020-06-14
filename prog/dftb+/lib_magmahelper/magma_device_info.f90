! helper routines for calls to the MAGMA library, Fortran layer on top of C routines
module device_info

  ! number of available GPUs on the system
  interface gpu_avail

    subroutine  gpu_avail4(max_ngpus) bind(c, name='MagmaGpuNum_Available')
      use iso_c_binding
      integer(c_int),intent(inout) :: max_ngpus
    end subroutine gpu_avail4

    subroutine  gpu_avail8(max_ngpus) bind(c, name='MagmaGpuNum_Available')
      use iso_c_binding
      integer(c_long),intent(inout) :: max_ngpus
    end subroutine gpu_avail8

  end interface gpu_avail

  ! number of GPUs requested by the code
  interface gpu_req

    subroutine  gpu_req4(req_ngpus) bind(c, name='MagmaGpuNum_Requested')
      use iso_c_binding
      integer(c_int),intent(inout) :: req_ngpus
    end subroutine gpu_req4

    subroutine  gpu_req8(req_ngpus) bind(c, name='MagmaGpuNum_Requested')
      use iso_c_binding
      integer(c_long),intent(inout) :: req_ngpus
    end subroutine gpu_req8

  end interface gpu_req

end module device_info
