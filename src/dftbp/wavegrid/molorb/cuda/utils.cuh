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

/**
 * @brief Multiplies a 3x3 matrix with a 3D vector.
 * @param mat    The 3x3 matrix.
 * @param vec    The 3D vector.
 * @param result The resulting 3D vector after multiplication.
 */
__device__ __forceinline__ void matmul3x3_vec(const double mat[3][3], const double vec[3], double result[3]) {
    for (int i = 0; i < 3; i++) {
        result[i] = mat[i][0] * vec[0] + mat[i][1] * vec[1] + mat[i][2] * vec[2];
    }
}

/**
 * @brief Multiplies a 3D vector with a 3x3 matrix.
 * @param vec    The 3D vector.
 * @param mat    The 3x3 matrix. 
 * @param result The resulting 3D vector after multiplication.
 */
__device__ __forceinline__ void vecmul3x3_mat(const double vec[3], const double mat[3][3], double result[3]) {
    for (int i = 0; i < 3; i++) {
        result[i] = vec[0] * mat[0][i] + vec[1] * mat[1][i] + vec[2] * mat[2][i];
    }
}


/**
 * @brief Folds coordinates into the unit cell by discarding non-fractional part in lattice vector multiples.
 * @param xyz         The 3D coordinates to be folded (input and output).
 * @param latVecs     The lattice vectors (3x3 matrix).
 * @param recVecs2p   Inverse of the lattice vecs (3x3 matrix, reciprocal lattice vectors divided by 2pi)
 */
__device__ __forceinline__ void foldCoordsIntoCell(
    double xyz[3], const double latVecs[3][3], const double recVecs2p[3][3]) {
    double frac[3];
    vecmul3x3_mat(xyz, recVecs2p, frac);

    for (int i = 0; i < 3; i++) {
        frac[i] -= floor(frac[i]);
    }

    matmul3x3_vec(latVecs, frac, xyz);
}
