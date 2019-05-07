#ifndef __TESTHELPERS_H__
#define __TESTHELPERS_H_

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Writes basic information for the autotest systsem
 *
 * \param n_atom[in] Number of atoms in the system.
 *
 * \param n_ext_charge[in]  Number of external charges in the system. Set it to zero,
 *     if no external charges are defined.
 *
 * \param mermin_energy[in]  Mermin free energy
 *
 * \param gradients[in]  Pointer to the array with the gradients. Pass NULL if not available.
 *
 * \param gross_charges[in]  Pointer to the array with the gross charges. Pass NULL if not
 *     available.
 *
 * \param ext_charges_gradients[in]  Pointer to the array with the gradients of the external charges
 *     Pass NULL if not available.
 */
void dftbp_write_autotest_tag(int n_atom, int n_ext_charge, double mermin_energy, double *gradients,
                              double *gross_charges, double * ext_charge_gradients);

#ifdef __cplusplus
}
#endif

#endif
