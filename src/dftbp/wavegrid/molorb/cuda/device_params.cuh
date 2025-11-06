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
    // Allocate nothing
    DeviceBuffer() = default;

    // Allocate only
    explicit DeviceBuffer(size_t count) { allocate(count); }

    // Allocate and copy from host
    DeviceBuffer(const T* host_ptr, size_t count) {
        allocate(count);
        copy_to_device(host_ptr, count);
    }

    // Destructor frees memory on device
    ~DeviceBuffer() { deallocate(); }

    // Disable copy constructor, assignment
    DeviceBuffer(const DeviceBuffer&)            = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    // Enable move semantics
    DeviceBuffer(DeviceBuffer&& other) noexcept :
        _devicePtr(other._devicePtr),
        _count(other._count) {
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
    void deallocate() {
        if (_devicePtr) {
            CHECK_CUDA(cudaFree(_devicePtr));
            _devicePtr = nullptr;
            _count     = 0;
        }
    }
    T*     _devicePtr = nullptr;
    size_t _count     = 0;
};

/*
 * We implement the lookup table as a 2D texture with a single float channel.
 * The radial Functions are stored row-wise.
 * We assume identical radial grids for all Orbitals.
 * GPUs have dedicated hardware for texture access & interpolation.
 */
class GpuLutTexture {
   public:
    GpuLutTexture(const double* lutData, int nPoints, int nOrbitals) {
        // Convert the Fortran passed doubles to floats
        size_t             totalValues = (size_t)nOrbitals * nPoints;
        std::vector<float> lutFloats(totalValues);
        for (size_t i = 0; i < totalValues; ++i)
            lutFloats[i] = static_cast<float>(lutData[i]);

        // Allocate memory on device
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        CHECK_CUDA(cudaMallocArray(&_lutArray, &channelDesc, nPoints, nOrbitals, 0));

        // Copy data to array.
        // We cannot use a 1d copy since our _lutArray is not a simple pointer
        // to contiguous global memory.
        CHECK_CUDA(cudaMemcpy2DToArray(_lutArray,  // dst array
            0,                        // no offset in dst
            0,                        // no offset in dst
            lutFloats.data(),         // src pointer
            nPoints * sizeof(float),  // src pitch (for alignment, bytes to next row)
            nPoints * sizeof(float),  // width in bytes
            nOrbitals,                // height (number of cached orbitals)
            cudaMemcpyHostToDevice));

        // Prepare texture object properties
        cudaResourceDesc resDesc{};
        resDesc.resType         = cudaResourceTypeArray;
        resDesc.res.array.array = _lutArray;

        cudaTextureDesc texDesc{};
        // OOB access clamped to edge values
        // (Should not occur if orb_cutoffs are set correctly)
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

    cudaTextureObject_t get() const { return _textureObject; }

    ~GpuLutTexture() { deallocate(); }

    // Disable copy constructor, assignment
    GpuLutTexture(const GpuLutTexture&)            = delete;
    GpuLutTexture& operator=(const GpuLutTexture&) = delete;

    // Enable move semantics
    GpuLutTexture(GpuLutTexture&& other) noexcept :
        _lutArray(other._lutArray), _textureObject(other._textureObject) {}

    GpuLutTexture& operator=(GpuLutTexture&& other) noexcept {
        if (this != &other) {
            deallocate();
            _lutArray = other._lutArray;
            _textureObject = other._textureObject;
            other._lutArray = nullptr;
            other._textureObject = 0;
        }

        return *this;
        
    }



   private:
    void deallocate() {
        if (_textureObject) cudaDestroyTextureObject(_textureObject);
        if (_lutArray) cudaFreeArray(_lutArray);
    }
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

    // Basis
    DeviceBuffer<int>    orb_angMoms;
    DeviceBuffer<double> orb_cutoffsSq;
    // Texture for radial lookup table
    GpuLutTexture orb_lut;

    // Eigenvectors
    DeviceBuffer<double>   eigVecsReal;
    DeviceBuffer<complexd> eigVecsCmpl;

    // Output (per-GPU batch buffer)
    DeviceBuffer<complexd> valueCmpl_out_batch;
    DeviceBuffer<double>   valueReal_out_batch;

    // Constructor handles all H2D allocation and copy
    DeviceData(const GridParams* grid, const SystemParams* system, const PeriodicParams* periodic,
        const StoBasisParams* basis, const CalculationParams* calc, int z_per_batch)
        : origin(grid->origin, 3),
          gridVecs(grid->gridVecs, 9),
          coords(system->coords, (size_t)3 * system->nAtom * system->nCell),
          species(system->species, system->nAtom),
          iStos(system->iStos, system->nSpecies + 1),
          orb_angMoms(basis->angMoms, basis->nOrbitals),
          orb_cutoffsSq(basis->cutoffsSq, basis->nOrbitals),
          orb_lut(basis->lutGridValues, basis->nLutPoints, basis->nOrbitals)
    {
        if (DEBUG) printf("Using radial lookup table with %d points for %d orbitals\n", basis->nLutPoints, basis->nOrbitals);

        if (calc->isRealInput) {
            eigVecsReal = DeviceBuffer<double>(calc->eigVecsReal, (size_t)system->nOrb * calc->nEigIn);
        } else {
            eigVecsCmpl = DeviceBuffer<complexd>(reinterpret_cast<const complexd*>(calc->eigVecsCmpl), (size_t)system->nOrb * calc->nEigIn);
            phases      = DeviceBuffer<complexd>(reinterpret_cast<const complexd*>(periodic->phases), (size_t)system->nCell * calc->nEigIn);
            kIndexes    = DeviceBuffer<int>(periodic->kIndexes, calc->nEigIn);
        }

        // Allocate the per-GPU batch output array
        size_t batch_buffer_size_elems = (size_t)grid->nPointsX * grid->nPointsY * z_per_batch * calc->nEigOut;
        if (calc->isRealOutput) {
            valueReal_out_batch = DeviceBuffer<double>(batch_buffer_size_elems);
        } else {
            valueCmpl_out_batch = DeviceBuffer<complexd>(batch_buffer_size_elems);
        }
    }
};

// Kernel parameters struct to simplify the argument list
struct DeviceKernelParams {
    // Grid
    const int     nPointsX, nPointsY;
    const double* origin;
    const double* gridVecs;

    // System
    int nAtom, nCell, nOrb;

    const double* coords;
    const int*    species;
    const int*    iStos;

    // Periodic boundary cond.
    const bool      isPeriodic;
    double          latVecs[3][3];
    double          recVecs2pi[3][3];
    const int*      kIndexes;
    const complexd* phases;

    // Basis
    const int nOrbitals;
    // Texture lookup tables
    cudaTextureObject_t lutTex;
    const double        inverseLutStep;
    const int*          orb_angMoms;
    const double*       orb_cutoffsSq;

    // Eigenvectors
    const int       nEig, nEig_per_pass;
    const double*   eigVecsReal;
    const complexd* eigVecsCmpl;

    // Output (batch pointers)
    double*   valueReal_out_batch;
    complexd* valueCmpl_out_batch;

    // batch offsets modified on each loop iteration
    int z_per_batch = 0;
    int z_offset_global = 0;

    // Constructor to initialize the parameters from host data
    // Batch-specific z-slices set to 0, need to be set in the loop before kernel launch.
    DeviceKernelParams(DeviceData& data, const GridParams* grid, const SystemParams* system,
        const PeriodicParams* periodic, const StoBasisParams* basis, const CalculationParams* calc,
        int nEig_per_pass_in) :
        // Grid
        origin(data.origin.get()),
        gridVecs(data.gridVecs.get()),
        nPointsX(grid->nPointsX),
        nPointsY(grid->nPointsY),
        // System
        nAtom(system->nAtom),
        nCell(system->nCell),
        nOrb(system->nOrb),
        coords(data.coords.get()),
        species(data.species.get()),
        iStos(data.iStos.get()),
        // Basis
        nOrbitals(basis->nOrbitals),
        orb_angMoms(data.orb_angMoms.get()),
        orb_cutoffsSq(data.orb_cutoffsSq.get()),
        lutTex(data.orb_lut.get()),
        inverseLutStep(basis->inverseLutStep),
        // Periodic boundary conditions
        isPeriodic(periodic->isPeriodic),
        kIndexes(data.kIndexes.get()),
        phases(data.phases.get()),
        // Eigenvectors
        nEig(calc->nEigIn),
        eigVecsReal(data.eigVecsReal.get()),
        eigVecsCmpl(data.eigVecsCmpl.get()),
        // Output batch buffers
        valueReal_out_batch(data.valueReal_out_batch.get()),
        valueCmpl_out_batch(data.valueCmpl_out_batch.get()),
        nEig_per_pass(nEig_per_pass_in)
    {
        // Periodic boundary conditions
        if (isPeriodic)
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    latVecs[i][j]    = periodic->latVecs[IDX2F(i, j, 3)];
                    recVecs2pi[i][j] = periodic->recVecs2pi[IDX2F(i, j, 3)];
            }
    }
};

