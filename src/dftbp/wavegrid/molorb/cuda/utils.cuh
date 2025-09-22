/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
#pragma once
#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

// Helper macro for robust CUDA calls
#define CHECK_CUDA(call)                                                         \
    do {                                                                         \
        cudaError_t __err = call;                                                \
        if (__err != cudaSuccess) {                                              \
            throw std::runtime_error(std::string("CUDA Error at ") +             \
                                     __FILE__ + ":" + std::to_string(__LINE__) + \
                                     " -> " + cudaGetErrorString(__err));        \
        }                                                                        \
    } while (0)

// Helper macros for column-major (Fortran-style) index calculations
// We cannot cast to explicit shape because dimensions need to be fixed at compile time
#define IDX2F(i, j, lda) ((j) * (size_t)(lda) + (i))
#define IDX3F(i, j, k, lda, ldb) (((k) * (size_t)(ldb) + (j)) * (size_t)(lda) + (i))
#define IDX4F(i, j, k, l, lda, ldb, ldc) ((((l) * (size_t)(ldc) + (k)) * (size_t)(ldb) + (j)) * (size_t)(lda) + (i))

// Enable additional print statements
#define DEBUG 0

