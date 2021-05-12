#include <magma.h>

// Function to count GPU

const int MAXGPUS = 117;  // large enough to fit all device IDs

// Obtain number of available of GPUs
void  dftbp_magma_get_gpus_available(int *max_ngpus)
{
    magma_device_t devices[MAXGPUS];
    magma_int_t number_of_gpus;

    magma_init();
    magma_getdevices(devices, MAXGPUS, &number_of_gpus);
    *max_ngpus = number_of_gpus;
    return;
}

// Number of GPUs requested by the code
void dftbp_magma_get_gpus_requested(int *req_ngpus)
{
    *req_ngpus= magma_num_gpus();
    return;
}
