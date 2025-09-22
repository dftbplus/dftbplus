/*-------------------------------------------------------------------------------------------------*
 *  DFTB+: general package for performing fast atomistic simulations                               *
 *  Copyright (C) 2006 - 2025  DFTB+ developers group                                              *
 *                                                                                                 *
 *  See the LICENSE file for terms of usage and distribution.                                      *
 *-------------------------------------------------------------------------------------------------*/
// Real spherical harmonics up to l=4.
// See also: spharmonics.F90
#pragma once

namespace RealTessYConsts {
    // l = 0
    constexpr double C0_0 = 0.28209479177387814; // 1/2 * sqrt(1/pi)

    // l = 1
    constexpr double C1 = 0.4886025119029199; // m =-1,0,1  1/2 * sqrt(3/pi)

    // l = 2
    constexpr double C2_0          = 0.31539156525252005; // m=0        1/4 * sqrt(5/pi)
    constexpr double C2_abs1_neg2  = 1.0925484305920792;  // m=-2,-1,1  1/2 * sqrt(15/pi)
    constexpr double C2_2          = 0.5462742152960396;  // m=2        1/4 * sqrt(15/pi)

    // l = 3
    constexpr double C3_0     = 0.3731763325901154;  // m=0     1/4 * sqrt(7/pi)
    constexpr double C3_abs1  = 0.4570457994644658;  // m=-1,1  1/4 * sqrt(21/(2*pi))
    constexpr double C3_neg2  = 2.890611442640554;   // m=-2    1/2 * sqrt(105/pi)
    constexpr double C3_2     = 1.445305721320277;   // m=2     1/4 * sqrt(105/pi)
    constexpr double C3_abs3  = 0.5900435899266435;  // m=-3,3  1/4 * sqrt(35/(2*pi))

    // l = 4
    constexpr double C4_0     = 0.10578554691520431; // m=0     3/16 * sqrt(1/pi)
    constexpr double C4_abs1  = 0.6690465435572892;  // m=-1,1  3/4 * sqrt(5/(2*pi))
    constexpr double C4_neg2  = 0.9461746957575601;  // m=-2    3/4 * sqrt(5/pi)
    constexpr double C4_2     = 0.47308734787878004; // m=2     3/8 * sqrt(5/pi)
    constexpr double C4_abs3  = 1.7701307697799304;  // m=-3,3  3/4 * sqrt(35/(2*pi))
    constexpr double C4_neg4  = 2.5033429417967046;  // m=-4    3/4 * sqrt(35/pi)
    constexpr double C4_4     = 0.6258357354491761;  // m=4     3/16 * sqrt(35/pi)

} // namespace RealTessYConsts


/**
 * @brief Computes real tesseral spherical harmonics Y_lm(r) up to l=4 without performing divisions.
 *
 * All Definitions in accordance with:
 * https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
 * See also (general set):
 * https://winter.group.shef.ac.uk/orbitron/atomic_orbitals/4f/4f_equations.html
 * 
 * 
 * @param l      Orbital quantum number
 * @param m      Magnetic quantum number
 * @param diff   Pointer to (x, y, z) vector
 * @param inv_r  Pre-calculated 1/r. (At origin, pass 0 to avoid NaNs)
 * @return The value of the real spherical harmonic.
 */
__device__ __forceinline__ double realTessY(int l, int m, const double* diff, double inv_r) {
    const double x = diff[0];
    const double y = diff[1];
    const double z = diff[2];

    const double x_r = x * inv_r;
    const double y_r = y * inv_r;
    const double z_r = z * inv_r;

    using namespace RealTessYConsts;

    switch (l) {
        case 0: // s orbital
            return C0_0;
        
        case 1: // p orbitals
            switch (m) {
                case -1: return C1 * y_r; // p_y
                case  0: return C1 * z_r; // p_z
                case  1: return C1 * x_r; // p_x
            }
            break;

        case 2: { // d orbitals
            const double x_r2 = x_r * x_r;
            const double y_r2 = y_r * y_r;
            const double z_r2 = z_r * z_r;
            switch (m) {
                case -2: return C2_abs1_neg2 * (x_r * y_r);   // d_xy
                case -1: return C2_abs1_neg2 * (y_r * z_r);   // d_yz
                case  0: return C2_0 * (3.0 * z_r2 - 1.0);    // d_z^2
                case  1: return C2_abs1_neg2 * (x_r * z_r);   // d_xz
                case  2: return C2_2 * (x_r2 - y_r2);         // d_x^2-y^2
            }
        } break;

        case 3: { // f orbitals
            const double x_r2 = x_r * x_r;
            const double y_r2 = y_r * y_r;
            const double z_r2 = z_r * z_r;
            switch (m) {
                case -3: return C3_abs3 * y_r * (3.0 * x_r2 - y_r2);   // f_y(3x^2-y^2)
                case -2: return C3_neg2 * (x_r * y_r * z_r);           // f_xyz
                case -1: return C3_abs1 * y_r * (5.0 * z_r2 - 1.0);    // f_y(5z^2-r^2) -> f_yz^2
                case  0: return C3_0 * z_r * (5.0 * z_r2 - 3.0);       // f_z(5z^2-3r^2) -> f_z^3
                case  1: return C3_abs1 * x_r * (5.0 * z_r2 - 1.0);    // f_x(5z^2-r^2) -> f_xz^2
                case  2: return C3_2 * z_r * (x_r2 - y_r2);            // f_z(x^2-y^2)
                case  3: return C3_abs3 * x_r * (x_r2 - 3.0 * y_r2);   // f_x(x^2-3y^2)
            }
        } break;

        case 4: { // g orbitals
            const double x_r2 = x_r * x_r;
            const double y_r2 = y_r * y_r;
            const double z_r2 = z_r * z_r;
            switch (m) {
                case -4: return C4_neg4 * (x_r * y_r * (x_r2 - y_r2));
                case -3: return C4_abs3 * (z_r * y_r * (3.0 * x_r2 - y_r2));
                case -2: return C4_neg2 * (x_r * y_r * (7.0 * z_r2 - 1.0));
                case -1: return C4_abs1 * (z_r * y_r * (7.0 * z_r2 - 3.0));
                case  0: return C4_0 * (35.0 * z_r2 * z_r2 - 30.0 * z_r2 + 3.0);
                case  1: return C4_abs1 * (z_r * x_r * (7.0 * z_r2 - 3.0));
                case  2: return C4_2 * ((x_r2 - y_r2) * (7.0 * z_r2 - 1.0));
                case  3: return C4_abs3 * (z_r * x_r * (x_r2 - 3.0 * y_r2));
                case  4: return C4_4 * (x_r2 * x_r2 - 6.0 * x_r2 * y_r2 + y_r2 * y_r2);
            }
        } break;
    }
    return 0.0;
}
