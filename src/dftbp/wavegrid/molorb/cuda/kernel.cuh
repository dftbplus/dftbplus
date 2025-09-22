/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
#ifndef KERNEL_CUH_
#define KERNEL_CUH_

#include <cuComplex.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Important:
 * These struct definitions must be mirrored at the c bound Fortran call site in offloaded.F90.
 */

// Calculation grid
typedef struct {
    const int     nPointsX;
    const int     nPointsY;
    const int     nPointsZ;
    const double* origin;    // [3]
    const double* gridVecs;  // [3][3]
} GridParams;

typedef struct {
    const int     nAtom;
    const int     nCell;
    const int     nSpecies;
    const int     nOrb;
    const double* coords;   // [3][nAtom][nCell]
    const int*    species;  // [nAtom]
    const int*    iStos;    // [nSpecies+1]
} SystemParams;

// Additional System information describing periodic boundary conditions
typedef struct {
    const bool             isPeriodic;
    const double*          latVecs;     // [3][3]
    const double*          recVecs2pi;  // [3][3]
    const int*             kIndexes;    // [nEig]
    const cuDoubleComplex* phases;      // [nCell][nEigIn]
} PeriodicParams;

// Basis parameters describing orbitals
typedef struct {
    const bool useRadialLut;
    const int  nStos;
    const int  nLutPoints;

    const double  inverseLutStep;
    const double* lutGridValues;  // [nStos][nLutPoints]

    const int     maxNPows;
    const int     maxNAlphas;
    const int*    angMoms;    // [nStos]
    const int*    nPows;      // [nStos]
    const int*    nAlphas;    // [nStos]
    const double* cutoffsSq;  // [nStos]
    const double* coeffs;     // [maxNPows][maxNAlphas][nStos]
    const double* alphas;     // [maxNAlphas][nStos]
} StoBasisParams;

// Coefficient Input and Control Flags
typedef struct {
    const bool             isRealInput;
    const bool             isRealOutput;
    const bool             calcAtomicDensity;
    const bool             calcTotalChrg;
    const int              nEigIn;
    const int              nEigOut;        // calcTotalChrg ? 1 : nEigIn
    const double*          eigVecsReal;    // [nOrb][nEigIn]
    const cuDoubleComplex* eigVecsCmpl;    // [nOrb][nEigIn]
    double*                valueReal_out;  // [nPointsX][nPointsY][nPointsZ][nEigOut]
    cuDoubleComplex*       valueCmpl_out;  // [nPointsX][nPointsY][nPointsZ][nEigOut]
} CalculationParams;

void evaluate_on_device_c(
        const GridParams* grid,
        const SystemParams* system,
        const PeriodicParams* periodic,
        const StoBasisParams* basis,
        const CalculationParams* calc);

#ifdef __cplusplus
}
#endif

#endif  // KERNEL_CUH_
