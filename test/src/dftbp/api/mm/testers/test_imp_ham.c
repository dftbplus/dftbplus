/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2024  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "dftbplus.h"
#include "testhelpers.h"

#define BASIS_SIZE 6
#define N_ATOMS 3

double dm[BASIS_SIZE][BASIS_SIZE];

/* hamiltonian to import: water molecule with 1e point charge  at (0, 0, 5)
 * point*/
const double hamiltonian[BASIS_SIZE][BASIS_SIZE] = {
    {-1.003919,  0.000000,  0.000000,  0.000000, -0.550573, -0.550573},
    { 0.000000, -0.457218,  0.000000,  0.000000, -0.309598,  0.309598},
    { 0.000000,  0.000000, -0.457218,  0.000000,  0.241885,  0.241885},
    { 0.000000,  0.000000,  0.000000, -0.457218,  0.000000,  0.000000},
    {-0.550573, -0.309598,  0.241885,  0.000000, -0.391745, -0.141157},
    {-0.550573,  0.309598,  0.241885,  0.000000, -0.141157, -0.391745}};

/*
Copy symmetric matrix with only uppert triangle stored.
*/
void copy_triu(const double *tril_src, double *dst)
{
  for (int i = 0; i < BASIS_SIZE; ++i)
  {
    for (int j = i; j < BASIS_SIZE; ++j)
    {
      dst[i * BASIS_SIZE + j] = tril_src[i * BASIS_SIZE + j];
      dst[j * BASIS_SIZE + i] = tril_src[i * BASIS_SIZE + j];
    }
  }
}

void dm_callback(void *aux_ptr, int iK, int iS, int *blacs_descr,
                 const void *blacs_data, ASI_matrix_descr_t *matrix_descr) {
  const double *dm_local = blacs_data;
  if (matrix_descr->storage_type!=ASI_STORAGE_TYPE_TRIL)
  {
    abort();
  }
  copy_triu(dm_local, &(dm[0][0]));
}

int set_h_callback(void *aux_ptr, int iK, int iS, int *blacs_descr,
                void *blacs_data, ASI_matrix_descr_t *matrix_descr) {

  if (matrix_descr->storage_type!=ASI_STORAGE_TYPE_TRIL)
  {
    abort();
  }

  double *h_local = blacs_data;
  copy_triu(&hamiltonian[0][0], h_local);

  return 1;
}

int main() {
  FILE *atf = fopen("autotest.tag", "w+");

  DftbPlus calculator;
  DftbPlusInput input;

  /* Coordinates in row major format, atomic units */
  double coords_h2o[] = {
    0.000000000000000000e+00,  0.000000000000000000e+00,  2.253725172195637505e-01,
    0.000000000000000000e+00,  1.442312678557651440e+00, -9.014881791521291987e-01,
    0.000000000000000000e+00, -1.442312678557651440e+00, -9.014881791521291987e-01,
  };

  int major, minor, patch;
  dftbp_api(&major, &minor, &patch);
  printf("API version %d.%d.%d\n", major, minor, patch);

  _Bool instsafe = dftbp_is_instance_safe();
  printf(instsafe ? "API is instance safe\n" : "API is NOT instance safe\n");

  /* Initialize DFTB+ input tree from input in external file */
  dftbp_init(&calculator, NULL);
  dftbp_get_input_from_file(&calculator, "dftb_in.h2o.hsd", &input);
  dftbp_process_input(&calculator, &input);

  int natom = dftbp_get_nr_atoms(&calculator);
  fprintf(atf, "natom       :integer:0:\n%d\n", natom);

  dftbp_set_coords(&calculator, coords_h2o);

  int basis_size = dftbp_get_basis_size(&calculator);
  fprintf(atf, "basis_size       :integer:0:\n%d\n", basis_size);

  _Bool is_hs_real = dftbp_is_hs_real(&calculator);
  assert(is_hs_real); // workaround for Intel(R) oneAPI 2024.0, as
		      // otherwise get value outside bool {0,1}
  fprintf(atf, "is_hs_real       :integer:0:\n%d\n", (is_hs_real ? 1 : 0));

  int n_spin = dftbp_get_nr_spin(&calculator);
  fprintf(atf, "n_spin       :integer:0:\n%d\n", n_spin);

  int n_kpts = dftbp_nr_kpoints(&calculator);
  fprintf(atf, "n_kpts       :integer:0:\n%d\n", n_kpts);

  dftbp_register_dm_callback(&calculator, dm_callback, 0);
  dftbp_register_set_h_callback(&calculator, set_h_callback, 0);

  assert(N_ATOMS == natom);
  assert(BASIS_SIZE == basis_size);
  assert(1 == n_kpts);
  assert(1 == n_spin);

  /* Evaluate energy */
  double mermin_energy;
  dftbp_get_energy(&calculator, &mermin_energy);
  fprintf(atf, "mermin_energy       :real:0:\n%f\n", mermin_energy);

  double EigSum = 0;
  for (int i = 0; i < BASIS_SIZE; ++i) {
    for (int j = 0; j < BASIS_SIZE; ++j) {
      EigSum += dm[i][j] * hamiltonian[i][j];
    }
  }
  fprintf(atf, "EigSum       :real:0:\n%f\n", EigSum);

  fclose(atf);
  return 0;
}
