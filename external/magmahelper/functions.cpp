#include <magma.h>

// count GPUs function

const int MAXGPUS = 117; // 117 is large enough to fit all device IDs

extern "C"  void  MagmaGpuNum_Available(int* max_ngpus)
{
    magma_device_t devices[MAXGPUS];
    magma_int_t number_of_gpus;
    magma_getdevices( devices, MAXGPUS, &number_of_gpus);
    *max_ngpus=number_of_gpus;
    return;
}

extern "C" void MagmaGpuNum_Requested(int* req_ngpus)
{
    *req_ngpus= magma_num_gpus();
    return;
}


