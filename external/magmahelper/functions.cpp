#include <magma.h>

// count GPUs function


extern "C"  void  MagmaGpuNum_Available(int* max_ngpus)
{
    magma_device_t devices[117]; // 117 is large enough to fit all device IDs
    magma_int_t number_of_gpus;
    magma_getdevices( devices, 117, &number_of_gpus);
    *max_ngpus=number_of_gpus;
    return;
}

extern "C" void MagmaGpuNum_Requested(int* req_ngpus)
{
    *req_ngpus= magma_num_gpus();
    return;
}


