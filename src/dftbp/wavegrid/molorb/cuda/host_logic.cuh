/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
#pragma once
#include "kernel.cuh"
#include "device_params.cuh"

#include <algorithm>
#include <cstdio>
#include <stdexcept>

// amount of shared memory set aside for nEig accumulators
constexpr float SHARED_MEM_FACTOR = 0.95f;
// max output array share of free global memory
// Leaves room for input arrays such as eigenvectors.
constexpr float GLOBAL_MEM_FACTOR = 0.80f;
// Should be a multiple of warp size 32
constexpr int THREADS_PER_BLOCK = 256;

// Debug kernel execution timing data
struct elapsedTime_ms {float kernel, d2h, everything;};

class GpuLaunchConfig {
   public:
    int    deviceId;            // CUDA device ID
    int    z_count;             // Number of Z-slices for this GPU
    int    z_start;             // Starting Z-slice for this GPU
    int    z_per_batch;         // Number of Z-slices to process per kernel launch
    int    nEig_per_pass;       // Number of eigenstates to accumulate in shared memory
    size_t shared_mem_for_pass; // Amount of shared memory per block for the accumulators
    // Kernel template parameters
    const bool isRealInput;
    const bool calcAtomicDensity;
    const bool calcTotalChrg;


    /** 
     * @brief Splits work among GPUs using Z-slices, and calculate gpu launch parameters based on available memory.
     *
     * @param deviceId      The CUDA GPU device ID.
     * @param numGpus       Total number of GPUs used for splitting the work.
     * @param grid          Pointer to the GridParams structure containing grid dimensions.
     * @param calc          Pointer to the CalculationParams structure containing calculation flags.
     * @return A GpuLaunchConfig structure with the calculated parameters for the specified GPU.
     */
    GpuLaunchConfig(int deviceId, int numGpus, const GridParams* grid, const CalculationParams* calc):
        deviceId(deviceId),
        isRealInput(calc->isRealInput),
        calcAtomicDensity(calc->calcAtomicDensity),
        calcTotalChrg(calc->calcTotalChrg)
    {
        // Evenly split Z-slices among GPUs
        distribute_z_slices(deviceId, numGpus, grid);
        // Query global memory and determine z_per_batch
        z_per_batch = determine_batch_size(grid, calc);
        // How many nEig can we fit in our shared shared memory accumulator?
        nEig_per_pass = determine_accumulator_size(grid, calc);
        // How much memory do we need for the accumulators?
        size_t accumulator_number_size = calc->isRealInput ? sizeof(double) : sizeof(complexd);
        shared_mem_for_pass = (size_t)nEig_per_pass * THREADS_PER_BLOCK * accumulator_number_size;
    }
    /**
     * @brief Prints a summary of the GPU configuration.
     * @param grid Pointer to the GridParams structure containing grid dimensions.
     * @param calc Pointer to the CalculationParams structure containing calculation parameters.
     */
    void print_summary(const GridParams* grid, const CalculationParams* calc) const {
        printf("\n--- GPU %d Configuration ---\n", deviceId);
        printf("  Z-slice workload: %d (from index %d to %d)\n", z_count, z_start,
            z_start + z_count - 1);
        printf("  Block size: %d threads, %zub shared mem per block, %d eigs of %d per pass\n", THREADS_PER_BLOCK,
            shared_mem_for_pass, nEig_per_pass, calc->nEigIn);
        size_t free_global_mem = get_available_global_bytes();
        size_t total_size_valueOut =
            (size_t)grid->nPointsX * grid->nPointsY * grid->nPointsZ * calc->nEigOut * sizeof(double);
        if (!calc->isRealOutput) total_size_valueOut *= 2;
        printf(" (Free device mem: %.2f GB, Grid size: %d x %d x %d (x %d eigs) = %.2f GB)\n", free_global_mem / 1e9,
            grid->nPointsX, grid->nPointsY, grid->nPointsZ, calc->nEigOut, total_size_valueOut / 1e9);
        printf("  Processing Z-slices in batches of %d\n", z_per_batch);
    }

   private:
    /**
     * @brief Determine the Z-slice range for this GPU based on device ID and total number of GPUs.
     * @param deviceId The CUDA GPU device ID (0 to numGpus-1).
     * @param numGpus Total number of GPUs used for splitting the work.
     * @param grid Pointer to the GridParams structure containing grid dimensions.
     */
    void distribute_z_slices(int deviceId, int numGpus, const GridParams* grid) {
        // Evenly distribute slices among GPUs
        int base_count = grid->nPointsZ / numGpus;
        int remainder = grid->nPointsZ % numGpus;
    
        // The first remainder GPUs (ID in 0...remainder-1) get additional slice
        z_count = base_count + (deviceId < remainder ? 1 : 0);
    
        // Calculate starting position
        z_start = deviceId * base_count + std::min(deviceId, remainder);
    }

    /**
     * @brief Determine the amount of available global memory on the device.
     * @return Available global memory in bytes, reduced by SHARED_MEM_FACTOR safety margin.
     */
    size_t get_available_global_bytes() const { 
        // Query available memory sizes with safety margins
        size_t         free_global_mem, total_global_mem;
        CHECK_CUDA(cudaMemGetInfo(&free_global_mem, &total_global_mem));
        size_t available_global = static_cast<size_t>(free_global_mem * GLOBAL_MEM_FACTOR);
        return available_global;
    }

    /**
     * @brief Determine the number of Z-slices that can be processed in a single batch,
     *        limited by available global memory for output arrays.
     * @param grid Pointer to the GridParams structure containing grid dimensions.
     * @param calc Pointer to the CalculationParams structure containing calculation flags.
     * @return Number of Z-slices that can be processed in one batch.
     */
    int determine_batch_size(const GridParams* grid, const CalculationParams* calc) const {
        size_t available_global   = get_available_global_bytes();
        size_t output_number_size = calc->isRealOutput ? sizeof(double) : sizeof(complexd);

        // Determine max Z-slices that can fit in available (global) memory
        size_t bytes_per_slice = (size_t)grid->nPointsX * grid->nPointsY * calc->nEigOut * output_number_size;
        int z_per_batch        = std::min(z_count, (int)(available_global / bytes_per_slice));
        if(!z_per_batch) throw std::runtime_error(
            "Insufficient global GPU memory available, unable to fit output array Z-slice.");
        return z_per_batch;
    }
    
    /**
     * @brief Determines the max amount of shared memory available to a block.
     * @return Available shared memory in bytes, reduced by SHARED_MEM_FACTOR safety margin.
     */
    size_t get_available_shared_bytes() const {
        // Query available memory sizes with safety margins
        cudaDeviceProp prop;
        CHECK_CUDA(cudaGetDeviceProperties(&prop, deviceId));
        size_t available_shared = prop.sharedMemPerBlock * SHARED_MEM_FACTOR;
        return available_shared;
    }

    /**
     * @brief Determine the number of eigenstates that each thread can accumulate in shared memory.
     * @param grid Pointer to the GridParams structure containing grid dimensions.
     * @param calc Pointer to the CalculationParams structure containing calculation flags.
     * @return Number of eigenstates that can be accumulated in shared memory per thread.
     */
    int determine_accumulator_size(const GridParams* grid, const CalculationParams* calc) const {
        size_t available_shared        = get_available_shared_bytes();
        size_t accumulator_number_size = calc->isRealInput ? sizeof(double) : sizeof(complexd);
 
        // Determine number of eigenstates that fit into shared memory
        int accumulator_size  = available_shared / (THREADS_PER_BLOCK * accumulator_number_size);
        accumulator_size      = std::min(calc->nEigIn, accumulator_size);
        if(!accumulator_size) throw std::runtime_error(
            "Insufficient shared GPU memory available, unable to fit eigenstate accumulators.");
        return accumulator_size;
    }


};





/**
 * @brief A simple GPU timer using CUDA events.
 *
 * Wraps CUDA event calls to provide a simple way to time GPU execution.
 * Use elapsed() to retrieve the accumulated time without stopping the timer.
 * The timer is initially stopped unless startNow=true is passed to the constructor.
 * Not thread safe.
 */
class GpuTimer {
public:

    /**
     * @brief Constructor. Creates the CUDA events, optionally starts the timer.
     * @param startNow If true, starts the timer immediately.
     */
    explicit GpuTimer(bool startNow = false) : _accumulated_ms(0.0f), _running(false) {
        CHECK_CUDA(cudaEventCreate(&_startEvent));
        CHECK_CUDA(cudaEventCreate(&_stopEvent));
        if (startNow) start();
        
    }

    ~GpuTimer() {
        cudaEventDestroy(_startEvent);
        cudaEventDestroy(_stopEvent);
    }
   
    /**
     * @brief Starts the timer. If already running, does nothing.
     */
    void start() {
        if (!_running) {
            CHECK_CUDA(cudaEventRecord(_startEvent));
            _running = true;
        }
    }
    /**
     * @brief Returns the total elapsed time, without stopping the timer.
     * If the timer is running, it records the current time and updates the accumulated time.
     * @return Elapsed time in milliseconds.
     */
    float elapsed_ms() {
        if (_running) {
            CHECK_CUDA(cudaEventRecord(_stopEvent));
            CHECK_CUDA(cudaEventSynchronize(_stopEvent));
            float elapsed;
            CHECK_CUDA(cudaEventElapsedTime(&elapsed, _startEvent, _stopEvent));
            CHECK_CUDA(cudaEventRecord(_startEvent));
            _accumulated_ms += elapsed;
        }
        return _accumulated_ms;
    }

    /**
     * @brief Stops the timer and returns the total elapsed time.
     * If the timer is not running, simply returns the accumulated time.
     * @return Total elapsed time in milliseconds.
     */
    float stop() {
        _accumulated_ms = elapsed_ms();
        _running = false;
        return _accumulated_ms;
    }

    // Disallow copy and assign
    GpuTimer(const GpuTimer&) = delete;
    GpuTimer& operator=(const GpuTimer&) = delete;

private:
    cudaEvent_t _startEvent{};
    cudaEvent_t _stopEvent{};
    float _accumulated_ms;
    bool _running;
};


// Function declarations
void copyD2H(void* d_src_ptr, void* h_dest_ptr, int nPointsX, int nPointsY, int nPointsZ, int z_per_batch,
    int z_offset_global, const CalculationParams* calc);
void dispatchKernel(const DeviceKernelParams* params, GpuLaunchConfig config, int grid_size);
elapsedTime_ms runBatchOnDevice(const GpuLaunchConfig& config, const GridParams* grid,
    const SystemParams* system, const PeriodicParams* periodic, const StoBasisParams* basis,
    const CalculationParams* calc);
