#include "utils.cuh"
// Maximum number of powers in STOs, for static array sizing
constexpr int STO_MAX_POWS = 16;

/**
 * @brief Computes the radial part of a Slater-type orbital (STO).
 * @param r         Distance from center
 * @param iL        Angular momentum
 * @param iOrb      Index of the orbital
 * @param nPows     Number of polynomial powers in the STO
 * @param nAlphas   Number of exponential terms in the STO
 * @param coeffs    Coefficients for the polynomial terms
 * @param alphas    Exponential decay constants
 * @param maxNPows  Maximum number of polynomial powers across all orbitals
 * @param maxNAlphas Maximum number of alphas across all orbitals
 * @return The radial value of the STO at distance r.
 */
__device__ __forceinline__ double getRadialValue(double r, int iL, int iOrb, int nPows, int nAlphas,
    const double* coeffs, const double* alphas, int maxNPows, int maxNAlphas) {
    double sto_tmp_pows[STO_MAX_POWS];

    double sto_tmp_rexp = 1.0;

    if (iL > 0 || r > 1.0e-12) {
        for (int p = 0; p < iL; ++p) {
            sto_tmp_rexp *= r;
        }
    }
    for (int ii = 0; ii < nPows; ++ii) {
        sto_tmp_pows[ii] = sto_tmp_rexp;
        sto_tmp_rexp *= r;
    }

    double radialVal = 0.0;
    for (int ii = 0; ii < nAlphas; ++ii) {
        double term = 0.0;
        for (int jj = 0; jj < nPows; ++jj) {
            term += coeffs[IDX3F(jj, ii, iOrb, maxNPows, maxNAlphas)] * sto_tmp_pows[jj];
        }
        radialVal += term * exp(alphas[IDX2F(ii, iOrb, maxNAlphas)] * r);
    }

    return radialVal;
}

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
