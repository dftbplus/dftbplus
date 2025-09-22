/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
#include <cuda_runtime.h>
#include <omp.h>
#include <thrust/complex.h>

#include <cstdio>
#include <stdexcept>

#include "kernel.cuh"
#include "host_logic.cuh"
#include "device_params.cuh"
#include "slater.cuh"
#include "spharmonics.cuh"
#include "utils.cuh"

// Avoid division by zero
constexpr double INV_R_EPSILON = 1.0e-12;

using complexd = thrust::complex<double>;

/**
 * @brief Main Molorb CUDA Kernel
 *
 * The kernel calculation closely follows the (easier to read) CPU implementation in parallel.F90.
 * Main differences arise due to the GPU architecture:
 * - We use shared memory to accumulate results per thread, which is much faster than global memory.
 * - Because shared memory is limited, we have to chunk the eigenstates into nEig_per_pass pieces.
 *   This means the spatial calculation is repeated for each chunk, but the accumulation is fast.
 * - The radial function lookup table uses hardware accelerated texture interpolation (if enabled).
 * - To avoid branching (dropped at compile time), we template the kernel 16 ways on boolean flags:
 *
 * @tparam isRealInput       whether real/complex eigenvectors are used (and adds phases).
 * @tparam useRadialLut      whether to use texture memory interpolation for the STO radial functions.
 * @tparam isPeriodic        enables folding of coords into the unit cell.
 * @tparam calcAtomicDensity squares the basis wavefunction contributions. In this case, the
 *                           occupation should be passed as the eigenvector.
 * @tparam calcTotalChrg     accumulates the density over all states in valueReal_out of shape (x,y,z,1).
 *                           Here, occupation should be passed by multiplying the eigenvecs with sqrt(occupation).
 *
 * @param p The DeviceKernelParams struct containing all necessary parameters for the kernel.
 */
template <bool isRealInput, bool calcAtomicDensity, bool calcTotalChrg, bool useRadialLut>
__global__ void evaluateKernel(const DeviceKernelParams p) {
    using AccumT = typename std::conditional<(isRealInput), double, complexd>::type;

    // Each thread gets its own private slice of the shared memory buffer for fast accumulation.
    // We have to chunk the eigenstates into nEig_per_pass due to size constraints.
    // (Cuda doesnt allow templating the shared memory type, so we simply recast it.)
    extern __shared__ char shared_workspace[];
    AccumT* point_results_pass = reinterpret_cast<AccumT*>(shared_workspace) + threadIdx.x * p.nEig_per_pass;

    // Map each thread to unique 1d index
    int idx_in_batch          = blockIdx.x * blockDim.x + threadIdx.x;
    int total_points_in_batch = p.nPointsX * p.nPointsY * p.z_per_batch;
    if (idx_in_batch >= total_points_in_batch) return;

    // Map 1d index to point in grid
    int i1        = idx_in_batch % p.nPointsX;
    int i2        = (idx_in_batch / p.nPointsX) % p.nPointsY;
    int i3_batch  = idx_in_batch / (p.nPointsX * p.nPointsY);
    int i3_global = i3_batch + p.z_offset_global;

    // Map point to global coordinates.
    double xyz[3];
    for (int i = 0; i < 3; ++i)
        xyz[i] = p.origin[i] + i1 * p.gridVecs[IDX2F(i, 0, 3)]
                      +        i2 * p.gridVecs[IDX2F(i, 1, 3)]
                      + i3_global * p.gridVecs[IDX2F(i, 2, 3)];

    // If periodic, fold into cell by discarding the non-fractional part in lattice vector multiples.
    if (p.isPeriodic) foldCoordsIntoCell(xyz, p.latVecs, p.recVecs2pi);

    double totChrgAcc = 0.0;
    // Loop over eigenstates in chunks that fit in shared memory
    for (int eig_base = 0; eig_base < p.nEig; eig_base += p.nEig_per_pass) {
        // Initialize the small, per-pass buffer for this thread
        for (int i = 0; i < p.nEig_per_pass; ++i) {
            point_results_pass[i] = AccumT(0.0);
        }

        // Since we run out of space in point_result_pass[], the spatial calculation
        // is repeated for each chunk of eigenstates.
        // This is to keep the accumulation in fast shared memory.
        for (int iCell = 0; iCell < p.nCell; ++iCell) {
            int orbital_idx_counter = 0;
            for (int iAtom = 0; iAtom < p.nAtom; ++iAtom) {
                int iSpecies = p.species[iAtom] - 1;

                double diff[3];
                for (int i = 0; i < 3; ++i) {
                    diff[i] = xyz[i] - p.coords[IDX3F(i, iAtom, iCell, 3, p.nAtom)];
                }
                double rr = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];

                for (int iOrb = p.iStos[iSpecies] - 1; iOrb < p.iStos[iSpecies + 1] - 1; ++iOrb) {
                    int iL = p.sto_angMoms[iOrb];
                    // Skip calculating all -l...+l orbitals if outside cutoff
                    if (rr > p.sto_cutoffsSq[iOrb]) {
                        orbital_idx_counter += 2 * iL + 1;
                        continue;
                    }
                    double r = sqrt(rr);

                    double radialVal;
                    if constexpr (useRadialLut) {
                        double lut_pos = 0.5 + r * p.inverseLutStep; // Add 0.5 to adress texel center (imagine pixels)
                        radialVal      = (double)tex2D<float>(p.lutTex, lut_pos, (float)iOrb + 0.5f);
                    } else {
                        radialVal = getRadialValue(r, iL, iOrb, p.sto_nPows[iOrb], p.sto_nAlphas[iOrb], p.sto_coeffs,
                            p.sto_alphas, p.maxNPows, p.maxNAlphas);
                    }

                    // precompute inverse used across several realTessY calls
                    double inv_r  = (r < INV_R_EPSILON) ? 0.0 : 1.0 / r;

                    for (int iM = -iL; iM <= iL; ++iM) {
                        double val = radialVal * realTessY(iL, iM, diff, inv_r);
                        if constexpr (calcAtomicDensity) val = val * val;

                        // Accumulate into the small shared memory buffer for the current chunk
                        for (int iEig_offset = 0; iEig_offset < p.nEig_per_pass; ++iEig_offset) {
                            int iEig = eig_base + iEig_offset;
                            if (iEig >= p.nEig) break;  // Don't go past the end on the last chunk
                            size_t eig_idx = IDX2F(orbital_idx_counter, iEig, p.nOrb);
                            if constexpr (isRealInput) {
                                point_results_pass[iEig_offset] += val * p.eigVecsReal[eig_idx];
                            } else {
                                point_results_pass[iEig_offset] +=
                                    val * p.phases[IDX2F(iCell, iEig, p.nCell)] * p.eigVecsCmpl[eig_idx];
                            }
                        }
                        orbital_idx_counter++;
                    }
                }
            }
        }

        // Write the complete nEig_per_pass chunk to global memory.
        for (int iEig_offset = 0; iEig_offset < p.nEig_per_pass; ++iEig_offset) {
            int iEig = eig_base + iEig_offset;
            if (iEig >= p.nEig) break;
            size_t out_idx = IDX4F(i1, i2, i3_batch, iEig, p.nPointsX, p.nPointsY, p.z_per_batch);
            if constexpr (isRealInput) {
                if constexpr (calcTotalChrg)
                    totChrgAcc += point_results_pass[iEig_offset] * point_results_pass[iEig_offset];
                else
                    p.valueReal_out_batch[out_idx] = point_results_pass[iEig_offset];

            } else {
                if constexpr (calcTotalChrg)
                    totChrgAcc += thrust::norm(point_results_pass[iEig_offset]);
                else
                    p.valueCmpl_out_batch[out_idx] = point_results_pass[iEig_offset];
            }
        }
    }

    // Density stored in first eig : (x,y,z,1)
    if constexpr (calcTotalChrg) {
        size_t out_idx = IDX4F(i1, i2, i3_batch, 0, p.nPointsX, p.nPointsY, p.z_per_batch);

        p.valueReal_out_batch[out_idx] = totChrgAcc;
    }
}


/**
 * @brief Dispatches the appropriate templated kernel based on input flags.
 *
 * Since the kernel is templated on the different calculation modes, we cannot simply pass booleans
 * at runtime and thus need this dispatch table to call the correct binary.
 *
 * @param params    Pointer to DeviceKernelParams struct containing all necessary parameters for the kernel.
 * @param config    The GpuLaunchConfig struct containing template flags and launch configuration.
 * @param grid_size Total number of thread blocks to launch.
 */
void dispatchKernel(const DeviceKernelParams* params, GpuLaunchConfig config, int grid_size) {
#define CALL_KERNEL(isReal, doAtomic, doChrg, useLut) \
    evaluateKernel<isReal, doAtomic, doChrg, useLut><<<grid_size, THREADS_PER_BLOCK, config.shared_mem_for_pass>>>(*params);

    int idx = (config.isRealInput       ? 1 : 0)
            + (config.calcAtomicDensity ? 2 : 0)
            + (config.calcTotalChrg     ? 4 : 0)
            + (config.useRadialLut      ? 8 : 0);

    // Refrain from compiling invalid combinations
    if(config.calcAtomicDensity && config.calcTotalChrg)
        throw std::runtime_error("Error: calcAtomicDensity and calcTotalChrg cannot both be true.\n");
    if(config.calcAtomicDensity && !config.isRealInput)
        throw std::runtime_error("Error: calcAtomicDensity requires real input vectors.\n");

    switch (idx) {
        case 0:  CALL_KERNEL(false, false, false, false); break;
        case 1:  CALL_KERNEL(true,  false, false, false); break;
        case 3:  CALL_KERNEL(true,  true,  false, false); break;
        case 4:  CALL_KERNEL(false, false, true,  false); break;
        case 5:  CALL_KERNEL(true,  false, true,  false); break;
        case 8:  CALL_KERNEL(false, false, false, true); break;
        case 9:  CALL_KERNEL(true,  false, false, true); break;
        case 11: CALL_KERNEL(true,  true,  false, true); break;
        case 12: CALL_KERNEL(false, false, true,  true); break;
        case 13: CALL_KERNEL(true,  false, true,  true); break;
        default: throw std::runtime_error("Error: invalid kernel configuration index.\n");
    }
#undef CALL_KERNEL
}



/**
 * @brief Entry point for evaluating molecular orbitals on a 3D grid using CUDA.
 *
 * C++ Host Interface (extern "C", then called from Fortran)
 * Launches the evaluation on available GPUs, splitting the work in Z-slices.
 *
 * @param grid     Pointer to GridParams struct containing grid parameters.
 * @param system   Pointer to SystemParams struct containing system parameters.
 * @param periodic Pointer to PeriodicParams struct containing periodic boundary parameters.
 * @param basis    Pointer to StoBasisParams struct containing STO basis parameters.
 * @param calc     Pointer to CalculationParams struct containing calculation parameters and output arrays.
 *                 Output is written to calc->valueReal_out or calc->valueCmpl_out depending on calc->isRealOutput.
 */
extern "C" void evaluate_on_device_c(const GridParams* grid, const SystemParams* system, const PeriodicParams* periodic,
    const StoBasisParams* basis, const CalculationParams* calc) {
    try {
        // We currently assume a hardcoded maximum for the number of powers.
        if (!basis->useRadialLut && basis->maxNPows > STO_MAX_POWS)
            throw std::runtime_error("Error: STO basis maxNPows exceeds STO_MAX_POWS.\n");
        if (calc->nEigIn * grid->nPointsX * grid->nPointsY * grid->nPointsZ == 0)
            throw std::runtime_error("Error: Zero-sized dimension in input parameters.\n");
        if (calc->calcTotalChrg && calc->nEigOut != 1)
            throw std::runtime_error("Error: When calculating total charge density, nEigOut must be 1.\n");
        if (!calc->calcTotalChrg && calc->nEigOut != calc->nEigIn)
            throw std::runtime_error("Error: nEigOut must match nEigIn unless calculating total charge density.\n");


        // Multi-GPU Setup
        int numGpus;
        CHECK_CUDA(cudaGetDeviceCount(&numGpus));
        if (numGpus == 0)
            throw std::runtime_error("No CUDA-enabled GPUs found. Unable to launch Kernel.\n");
        printf("Wavegrid: Found %d CUDA-enabled GPUs.\n", numGpus);

#ifndef _OPENMP
        if (numGpus > 1) {
            fprintf(stderr, "\nWARNING: Code not compiled with OpenMP support (-fopenmp). Falling back to single-GPU mode.\n");
            numGpus = 1;
            printf("Running on GPU 0 only.\n");
        }
#endif
        elapsedTime_ms timings, threadTimings;
        // Use OMP to split across available GPUs
        // This works irrespective of the number of threads set in OMP_NUM_THREADS.
        #pragma omp parallel num_threads(numGpus)
        {
            int deviceId = omp_get_thread_num();
            GpuLaunchConfig config(deviceId, numGpus, grid, calc, basis->useRadialLut);
            if(DEBUG) config.print_summary(grid, calc);

            threadTimings = runBatchOnDevice(config, grid, system, periodic, basis, calc);
            if (deviceId == 0) timings = threadTimings;
        } // End of omp parallel region

        if (DEBUG) {
            printf("\nGPU 0 execution time: %.1f ms\n", timings.everything);
            float kernel_share = timings.kernel / timings.everything;
            printf(" -> Kernel:   %.1f ms (%.1f%%)\n", timings.kernel, kernel_share * 100.0f);
            printf(" -> D2H copy: %.1f ms (%.1f%%)\n", timings.d2h, (1.0f - kernel_share) * 100.0f);
        }

    // Since the Fortran library does not currently support exceptions, we simply
    // terminate the program upon encountering any.
    } catch (const std::exception& e) {
        fprintf(stderr, "Error during CUDA offloading, in evaluate_on_device_c:");
        fprintf(stderr, "%s\n", e.what());
        exit(EXIT_FAILURE);
    }
}



