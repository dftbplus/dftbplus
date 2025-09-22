/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
// This file contains the kernel parameter struct and
// helpers to manage allocating and copying data to the GPU.
#pragma once

#include <memory>
#include <stdexcept>
#include <vector>

#include <cuda_runtime.h>
#include <thrust/complex.h>

#include "kernel.cuh"
#include "utils.cuh"


using complexd = thrust::complex<double>;


/*
 * A simple RAII wrapper for device memory.
 * Use get() to retrieve the raw pointer.
 */
template <typename T>
class DeviceBuffer {
   public:
    DeviceBuffer() = default;

    explicit DeviceBuffer(size_t count) { allocate(count); }

    // Allocate and copy from host
    DeviceBuffer(const T* host_ptr, size_t count) {
        allocate(count);
        copy_to_device(host_ptr, count);
    }

    // Destructor automatically frees the memory
    ~DeviceBuffer() { deallocate(); }

    void deallocate() {
        if (_devicePtr) {
            CHECK_CUDA(cudaFree(_devicePtr));
            _devicePtr = nullptr;
            _count     = 0;
        }
    }

    // Assign a new size and copy from host
    void assign(const T* host_ptr, size_t count) {
        deallocate();
        allocate(count);
        copy_to_device(host_ptr, count);
    }

    // Disable copy semantics
    DeviceBuffer(const DeviceBuffer&)            = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    // Enable move semantics
    DeviceBuffer(DeviceBuffer&& other) noexcept : _devicePtr(other._devicePtr), _count(other._count) {
        other._devicePtr = nullptr;
        other._count     = 0;
    }
    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept {
        if (this != &other) {
            deallocate();
            _devicePtr       = other._devicePtr;
            _count           = other._count;
            other._devicePtr = nullptr;
            other._count     = 0;
        }
        return *this;
    }

    T*       get() { return _devicePtr; }
    const T* get() const { return _devicePtr; }
    size_t   size() const { return _count; }

    void copy_to_host(T* host_ptr, size_t count_to_copy) const {
        if (count_to_copy > _count) 
            throw std::runtime_error("Error: trying to copy more elements than buffer contains.");

        CHECK_CUDA(cudaMemcpy(host_ptr, _devicePtr, count_to_copy * sizeof(T), cudaMemcpyDeviceToHost));
    }

    void copy_to_device(const T* host_ptr, size_t count_to_copy) {
        if (!_devicePtr) 
            throw std::runtime_error("Error: device pointer is null. Cannot copy to device.");
        if (count_to_copy > _count)
            throw std::runtime_error("Error: trying to copy more elements than buffer contains.");

        CHECK_CUDA(cudaMemcpy(_devicePtr, host_ptr, count_to_copy * sizeof(T), cudaMemcpyHostToDevice));
    }

   private:
    void allocate(size_t count) {
        _count = count;
        if (count > 0) {
            CHECK_CUDA(cudaMalloc(&_devicePtr, count * sizeof(T)));
        } else {
            _devicePtr = nullptr;
        }
    }
    T*     _devicePtr = nullptr;
    size_t _count     = 0;
};

/*
 * We implement the LUT as a 2D texture with a single float channel.
 * The radial Functions are stored row-wise.
 * We assume identical radial grids for all STOs.
 * GPUs have dedicated hardware for texture access & interpolation.
 */
class GpuLutTexture {
   public:
    GpuLutTexture(const double* lutData, int nPoints, int nStos) {
        // Convert the Fortran passed doubles to floats
        size_t             totalValues = (size_t)nStos * nPoints;
        std::vector<float> lutFloats(totalValues);
        for (size_t i = 0; i < totalValues; ++i)
            lutFloats[i] = static_cast<float>(lutData[i]);

        // Allocate memory on device
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        CHECK_CUDA(cudaMallocArray(&_lutArray, &channelDesc, nPoints, nStos, 0));

        // Copy data to array
        CHECK_CUDA(cudaMemcpy2DToArray(_lutArray,  // dst array
            0,                        // no offset in dst
            0,                        // no offset in dst
            lutFloats.data(),         // src pointer
            nPoints * sizeof(float),  // src pitch (for alignment, bytes to next row)
            nPoints * sizeof(float),  // width in bytes
            nStos,                    // height (number of cached stos)
            cudaMemcpyHostToDevice));

        // Prepare texture object properties
        cudaResourceDesc resDesc{};
        resDesc.resType         = cudaResourceTypeArray;
        resDesc.res.array.array = _lutArray;

        cudaTextureDesc texDesc{};
        // OOB access clamped to edge values
        // (Should not occur if sto_cutoffs are set correctly)
        texDesc.addressMode[0] = cudaAddressModeClamp;
        texDesc.addressMode[1] = cudaAddressModeClamp;
        // Enable linear interpolation
        texDesc.filterMode = cudaFilterModeLinear;
        // Do not normalize the lut values
        texDesc.readMode = cudaReadModeElementType;
        // Access using texel coords, i.e. [0, N-1]
        // (Imagine pixels, add 0.5f to get the center.)
        texDesc.normalizedCoords = 0;

        // Create texture object
        CHECK_CUDA(cudaCreateTextureObject(&_textureObject, &resDesc, &texDesc, nullptr));
    }

    ~GpuLutTexture() {
        if (_textureObject) cudaDestroyTextureObject(_textureObject);
        if (_lutArray) cudaFreeArray(_lutArray);
    }

    // Disable copy
    GpuLutTexture(const GpuLutTexture&)            = delete;
    GpuLutTexture& operator=(const GpuLutTexture&) = delete;

    cudaTextureObject_t get() const { return _textureObject; }

   private:
    cudaArray_t         _lutArray      = nullptr;
    cudaTextureObject_t _textureObject = 0;
};

// Manages Gpu memory allocation and H2D copy
struct DeviceData {
    // Grid
    DeviceBuffer<double> origin;
    DeviceBuffer<double> gridVecs;

    // System
    DeviceBuffer<double> coords;
    DeviceBuffer<int>    species;
    DeviceBuffer<int>    iStos;

    // Periodic
    DeviceBuffer<int>      kIndexes;
    DeviceBuffer<complexd> phases;

    // STO Basis
    DeviceBuffer<int>    sto_angMoms;
    DeviceBuffer<int>    sto_nPows;
    DeviceBuffer<int>    sto_nAlphas;
    DeviceBuffer<double> sto_cutoffsSq;
    DeviceBuffer<double> sto_coeffs;
    DeviceBuffer<double> sto_alphas;
    // Texture for radial LUT
    std::unique_ptr<GpuLutTexture> sto_lut;

    // Eigenvectors
    DeviceBuffer<double>   eigVecsReal;
    DeviceBuffer<complexd> eigVecsCmpl;

    // Output (per-GPU batch buffer)
    DeviceBuffer<complexd> d_valueCmpl_out_batch;
    DeviceBuffer<double>   d_valueReal_out_batch;

    // Constructor handles all H2D allocation and copy
    DeviceData(const GridParams* grid, const SystemParams* system, const PeriodicParams* periodic,
        const StoBasisParams* basis, const CalculationParams* calc, int z_per_batch)
        : origin(grid->origin, 3),
          gridVecs(grid->gridVecs, 9),
          coords(system->coords, (size_t)3 * system->nAtom * system->nCell),
          species(system->species, system->nAtom),
          iStos(system->iStos, system->nSpecies + 1),
          sto_angMoms(basis->angMoms, basis->nStos),
          sto_cutoffsSq(basis->cutoffsSq, basis->nStos) {
        if (basis->useRadialLut) {
            if (DEBUG) printf("Using radial LUT with %d points for %d STOs\n", basis->nLutPoints, basis->nStos);
            sto_lut = std::unique_ptr<GpuLutTexture>(
                new GpuLutTexture(basis->lutGridValues, basis->nLutPoints, basis->nStos));
        } else {
            if (DEBUG) printf("Using direct STO evaluation for %d STOs\n", basis->nStos);
            sto_nPows.assign(basis->nPows, basis->nStos);
            sto_nAlphas.assign(basis->nAlphas, basis->nStos);
            sto_coeffs.assign(basis->coeffs, (size_t)basis->maxNPows * basis->maxNAlphas * basis->nStos);
            sto_alphas.assign(basis->alphas, (size_t)basis->maxNAlphas * basis->nStos);
        }

        if (calc->isRealInput) {
            eigVecsReal.assign(calc->eigVecsReal, (size_t)system->nOrb * calc->nEigIn);
        } else {
            eigVecsCmpl.assign(
                reinterpret_cast<const complexd*>(calc->eigVecsCmpl), (size_t)system->nOrb * calc->nEigIn);
            phases.assign(reinterpret_cast<const complexd*>(periodic->phases), (size_t)system->nCell * calc->nEigIn);
            kIndexes.assign(periodic->kIndexes, calc->nEigIn);
        }

        // Per-GPU batch buffer for the output
        size_t batch_buffer_size_elems = (size_t)grid->nPointsX * grid->nPointsY * z_per_batch * calc->nEigOut;
        if (calc->isRealOutput) {
            d_valueReal_out_batch = DeviceBuffer<double>(batch_buffer_size_elems);
        } else {
            d_valueCmpl_out_batch = DeviceBuffer<complexd>(batch_buffer_size_elems);
        }
    }
};

// Kernel parameters struct to simplify the argument list
struct DeviceKernelParams {
    // Grid
    int           nPointsX, nPointsY, z_per_batch, z_offset_global;
    const double* origin;
    const double* gridVecs;

    // System
    int nAtom, nCell, nOrb;

    const double* coords;
    const int*    species;
    const int*    iStos;

    // Periodic boundary cond.
    bool            isPeriodic;
    double          latVecs[3][3];
    double          recVecs2pi[3][3];
    const int*      kIndexes;
    const complexd* phases;

    // STO Basis
    int nStos, maxNPows, maxNAlphas;
    // Texture LUTs
    cudaTextureObject_t lutTex;
    double              inverseLutStep;
    // STO parameters
    const int*    sto_angMoms;
    const int*    sto_nPows;
    const int*    sto_nAlphas;
    const double* sto_cutoffsSq;
    const double* sto_coeffs;
    const double* sto_alphas;

    // Eigenvectors
    int             nEig, nEig_per_pass;
    const double*   eigVecsReal;
    const complexd* eigVecsCmpl;

    // Output (batch pointers)
    double*   valueReal_out_batch;
    complexd* valueCmpl_out_batch;

    // Constructor to initialize the parameters from host data
    // Batch-specific parameters are initialized to zero or nullptr,
    // and need to be set in the loop before kernel launch.
    DeviceKernelParams(DeviceData& data, const GridParams* grid, const SystemParams* system,
        const PeriodicParams* periodic, const StoBasisParams* basis, const CalculationParams* calc,
        int nEig_per_pass_in) {
        // Grid
        origin   = data.origin.get();
        gridVecs = data.gridVecs.get();
        nPointsX = grid->nPointsX;
        nPointsY = grid->nPointsY;

        // System
        nAtom   = system->nAtom;
        nCell   = system->nCell;
        nOrb    = system->nOrb;
        coords  = data.coords.get();
        species = data.species.get();
        iStos   = data.iStos.get();

        // STO Basis
        nStos         = basis->nStos;
        sto_angMoms   = data.sto_angMoms.get();
        sto_cutoffsSq = data.sto_cutoffsSq.get();

        if (basis->useRadialLut) {
            lutTex         = data.sto_lut->get();
            inverseLutStep = basis->inverseLutStep;
        } else {
            maxNPows    = basis->maxNPows;
            maxNAlphas  = basis->maxNAlphas;
            sto_nPows   = data.sto_nPows.get();
            sto_nAlphas = data.sto_nAlphas.get();
            sto_coeffs  = data.sto_coeffs.get();
            sto_alphas  = data.sto_alphas.get();
        }

        // Periodic boundary conditions
        isPeriodic = periodic->isPeriodic;
        kIndexes   = data.kIndexes.get();
        phases     = data.phases.get();
        if (isPeriodic)
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    latVecs[i][j]    = periodic->latVecs[IDX2F(i, j, 3)];
                    recVecs2pi[i][j] = periodic->recVecs2pi[IDX2F(i, j, 3)];
            }

        // Eigenvectors
        nEig        = calc->nEigIn;
        eigVecsReal = data.eigVecsReal.get();
        eigVecsCmpl = data.eigVecsCmpl.get();

        // Output batch buffers
        valueReal_out_batch = data.d_valueReal_out_batch.get();
        valueCmpl_out_batch = data.d_valueCmpl_out_batch.get();

        nEig_per_pass = nEig_per_pass_in;
        // Batch-specific kernel config to be updated in the loop
        z_per_batch     = 0;
        z_offset_global = 0;
    }
};

