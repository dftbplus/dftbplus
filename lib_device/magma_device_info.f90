module device_info
    interface gpu_avail
        subroutine  gpu_avail4(max_ngpus) bind(c, name='MagmaGpuNum_Available')
            use iso_c_binding
            integer(c_int),intent(inout)     ::  max_ngpus
        end subroutine

        subroutine  gpu_avail8(max_ngpus) bind(c, name='MagmaGpuNum_Available')
            use iso_c_binding
            integer(c_long),intent(inout)    ::  max_ngpus
        end subroutine
    end interface gpu_avail
    
    interface gpu_req
        subroutine  gpu_req4(req_ngpus) bind(c, name='MagmaGpuNum_Requested')
            use iso_c_binding
            integer(c_int),intent(inout)     ::  req_ngpus
        end subroutine

        subroutine  gpu_req8(req_ngpus) bind(c, name='MagmaGpuNum_Requested')
            use iso_c_binding
            integer(c_long),intent(inout)    ::  req_ngpus
        end subroutine
    end interface gpu_req

end module device_info
