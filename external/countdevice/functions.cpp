#include <magma.h>

// count GPUs function


extern "C"  void  MagmaGpuNum_Available(int* max_ngpus)
{
    int max_num_gpus=200; //200 large enough to fit all device IDs
    magma_device_t devices[max_num_gpus];
    magma_int_t number_of_gpus;
    magma_getdevices( devices, max_num_gpus, &number_of_gpus);
    *max_ngpus=number_of_gpus;
    return;
}

extern "C" void MagmaGpuNum_Requested(int* req_ngpus)
{
    *req_ngpus= magma_num_gpus();
    return;
}


