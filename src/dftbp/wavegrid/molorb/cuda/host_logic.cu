/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
// The functions contained herein handle the high-level flow of data transfer,
// batchwise kernel execution and d2h-copy on a single GPU.
#include "host_logic.cuh"
#include <cuda_runtime.h>


/**
 * @brief Evaluates the assigned batch on a single GPU.
 *
 * Copies Data to the device, launches the kernel in a loop over Z-slices,
 * and copies the results back to the host.
 * Also returns elapsed time during kernel execution and D2H copy.
 *
 * @param config   The GpuLaunchConfig struct containing launch configuration and template flags.
 * @param grid     Pointer to GridParams struct containing grid parameters.
 * @param system   Pointer to SystemParams struct containing system parameters.
 * @param periodic Pointer to PeriodicParams struct containing periodic boundary parameters.
 * @param basis    Pointer to StoBasisParams struct containing STO basis parameters.
 * @param calc     Pointer to CalculationParams struct containing calculation parameters and output arrays.
 * @return An elapsedTime_ms struct containing the time spent in kernel execution and D2H copy.
 */
elapsedTime_ms runBatchOnDevice(const GpuLaunchConfig& config, const GridParams* grid,
    const SystemParams* system, const PeriodicParams* periodic, const StoBasisParams* basis,
    const CalculationParams* calc) {
    if (!config.z_count) return {0.0f, 0.0f, 0.0f};

    CHECK_CUDA(cudaSetDevice(config.deviceId));

    GpuTimer everything_timer(true), kernel_timer, d2h_timer;

    // Device allocation and H2D transfer
    DeviceData device_data(grid, system, periodic, basis, calc, config.z_per_batch);

    // Populate Kernel Parameter struct
    DeviceKernelParams deviceParams(device_data, grid, system, periodic, basis, calc, config.nEig_per_pass);

    // Process task in in batches (Z-slices)
    for (int z_offset = 0; z_offset < config.z_count; z_offset += config.z_per_batch) {
        deviceParams.z_per_batch = std::min(config.z_per_batch, config.z_count - z_offset);
        deviceParams.z_offset_global = config.z_start + z_offset;  // required to calculate coordinates in kernel

        int total_points_in_batch = grid->nPointsX * grid->nPointsY * deviceParams.z_per_batch;
        int grid_size             = (total_points_in_batch + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

        kernel_timer.start();
        dispatchKernel(&deviceParams, config, grid_size);
        kernel_timer.stop();

        d2h_timer.start();
        copyD2H(calc->isRealOutput ? (void*)device_data.valueReal_out_batch.get()
                                   : (void*)device_data.valueCmpl_out_batch.get(),
            calc->isRealOutput ? (void*)calc->valueReal_out : (void*)calc->valueCmpl_out, grid->nPointsX,
            grid->nPointsY, grid->nPointsZ, deviceParams.z_per_batch, deviceParams.z_offset_global, calc);
        d2h_timer.stop();
    }

    // Ensure we finished before returning
    CHECK_CUDA(cudaDeviceSynchronize());
    return {kernel_timer.elapsed_ms(), d2h_timer.elapsed_ms(), everything_timer.elapsed_ms()};
}

/** 
 * @brief Copies a batch of computed eigenstates from device to host memory.
 *
 * The D2H copy will automatically block/ synchronize the kernel for this batch.
 * This could be improved by using streams / cudaMemcpyAsync, but currently is not a bottleneck.
 * The output array is of fortran shape (x,y,z,nEigOut), thus z-slices are not contiguous.
 * We slice on Z instead of nEigOut to save on a little computation in the kernel.
 *
 * @param d_src_ptr        Pointer to the source data on the device (GPU).
 * @param h_dest_ptr       Pointer to the destination data on the host (CPU).
 * @param nPointsX         Number of grid points in the X dimension.
 * @param nPointsY         Number of grid points in the Y dimension.
 * @param nPointsZ         Total number of grid points in the Z dimension.
 * @param z_per_batch      Number of Z-slices processed in the current batch.
 * @param z_offset_global  Global Z-offset indicating where this batch fits in the full grid.
 * @param calc             Pointer to the CalculationParams structure containing calculation flags.
 */
void copyD2H(void* d_src_ptr, void* h_dest_ptr, int nPointsX, int nPointsY, int nPointsZ, int z_per_batch,
    int z_offset_global, const CalculationParams* calc) {
    size_t output_num_size = calc->isRealOutput ? sizeof(double) : sizeof(complexd);
    size_t host_plane_size    = (size_t)nPointsZ * nPointsY * nPointsX * output_num_size;
    size_t device_plane_size  = (size_t)z_per_batch * nPointsY * nPointsX * output_num_size;

    // Strided Memcpy3D could be used to squash this loop.
    for (int iEig = 0; iEig < calc->nEigOut; ++iEig) {
        // From: iEig-th slice of GPU batch buffer
        ptrdiff_t d_offset_bytes = (ptrdiff_t)(iEig * device_plane_size);

        // To: Global Z-position in the iEig-th slice of host buffer
        ptrdiff_t h_offset_bytes = (ptrdiff_t)(iEig * host_plane_size + ((size_t)z_offset_global * nPointsY * nPointsX) * output_num_size);

        CHECK_CUDA(cudaMemcpy((char*)h_dest_ptr + h_offset_bytes, (char*)d_src_ptr + d_offset_bytes, device_plane_size,
            cudaMemcpyDeviceToHost));
    }
}


