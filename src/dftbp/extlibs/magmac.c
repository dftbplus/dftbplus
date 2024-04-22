/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2023  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

/**
 *  Functions to count or request GPUs
 */

#include <magma.h>

/**
 *  Large enough to fit all device IDs
 */
const int MAXGPUS = 117;

/**
 *   Obtain number of available GPUs
 */
void  dftbp_extlibs_magma_get_gpus_available(int *max_ngpus)
{
    magma_device_t devices[MAXGPUS];
    magma_int_t number_of_gpus;

    magma_getdevices(devices, MAXGPUS, &number_of_gpus);
    *max_ngpus = number_of_gpus;
    return;
}

/**
 * Number of GPUs requested by the code
 */
void dftbp_extlibs_magma_get_gpus_requested(int *req_ngpus)
{
    *req_ngpus= magma_num_gpus();
    return;
}
