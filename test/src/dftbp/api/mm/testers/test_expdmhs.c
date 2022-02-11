/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2022  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>
#include "dftbplus.h"
#include "testhelpers.h"

#define BASIS_SIZE 8
#define N_KPTS 4
#define N_SPIN 1

double complex dm[N_KPTS][N_SPIN][BASIS_SIZE][BASIS_SIZE];
double complex overlap[BASIS_SIZE][BASIS_SIZE];
double complex hamiltonian[BASIS_SIZE][BASIS_SIZE];

void dm_callback(void *aux_ptr, int iK, int iS, int *blacs_descr, void *blacs_data)
{
  double complex *dm_local = blacs_data;
  for (int i = 0; i < BASIS_SIZE; ++i)
  for (int j = 0; j < BASIS_SIZE; ++j)
    dm[iK-1][iS-1][i][j] = dm_local[i*BASIS_SIZE + j];
}

void s_callback(void *aux_ptr, int *blacs_descr, void *blacs_data)
{
  double complex *s_local = blacs_data;
  for (int i = 0; i < BASIS_SIZE; ++i)
  for (int j = 0; j < BASIS_SIZE; ++j)
    overlap[i][j] = s_local[i*BASIS_SIZE + j];
}

void h_callback(void *aux_ptr, int *blacs_descr, void *blacs_data)
{
  double complex *h_local = blacs_data;
  for (int i = 0; i < BASIS_SIZE; ++i)
  for (int j = 0; j < BASIS_SIZE; ++j)
    hamiltonian[i][j] = h_local[i*BASIS_SIZE + j];
}

void print_matrix(FILE *f, const double complex *m)
{
    for (int i = 0; i < BASIS_SIZE; ++i)
    {
      for (int j = 0; j < BASIS_SIZE; ++j)
      {
        const double complex dc = m[i*BASIS_SIZE + j];
        fprintf(f, "%f %f\n", creal(dc), cimag(dc));
      }
    }
}

int main()
{
  FILE *atf = fopen("autotest.tag", "w+");
  fprintf(atf, "dummy_var       :real:0:\n%f\n",123.0);
  
  DftbPlus calculator;
  DftbPlusInput input;

  /* Coordinates in row major format, atomic units */
  double coords_si2[] = {
    0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
    2.2639291987021915, 2.4639291987021915, 2.5639291987021915
  };

  /* Lattice vectors in row major format, atomic units */
  double latvecs_si2[] = {
    5.2278583974043830, 5.1278583974043830, 0.0000000000000000,
    0.0000000000000000, 5.3278583974043830, 5.1278583974043830,
    5.1278583974043830, 0.0000000000000000, 5.4278583974043830
  };

  int major, minor, patch;
  dftbp_api(&major, &minor, &patch);
  printf("API version %d.%d.%d\n", major, minor, patch);

  _Bool instsafe = dftbp_is_instance_safe();
  printf(instsafe ? "API is instance safe\n" : "API is NOT instance safe\n");

  /* Initialize DFTB+ input tree from input in external file */
  dftbp_init(&calculator, NULL);
  dftbp_get_input_from_file(&calculator, "dftb_in.Si2.hsd", &input);
  dftbp_process_input(&calculator, &input);
  
  int natom = dftbp_get_nr_atoms(&calculator);
  fprintf(atf, "natom       :integer:0:\n%d\n",natom);
  
  dftbp_set_coords_and_lattice_vecs(&calculator, coords_si2, latvecs_si2);
  
  int basis_size = dftbp_get_basis_size(&calculator);
  fprintf(atf, "basis_size       :integer:0:\n%d\n", basis_size);
  
  _Bool is_hs_real = dftbp_is_hs_real(&calculator);
  fprintf(atf, "is_hs_real       :integer:0:\n%d\n", (is_hs_real ? 1 : 0) );
  
  int n_spin = dftbp_get_nr_spin(&calculator);
  fprintf(atf, "n_spin       :integer:0:\n%d\n", n_spin);
  
  int n_kpts = dftbp_nr_kpoints(&calculator);
  fprintf(atf, "n_kpts       :integer:0:\n%d\n", n_kpts);
  
  double kweights[N_KPTS];
  dftbp_get_kweights(&calculator, kweights);
  fprintf(atf, "kweights       :real:1:4\n");
  fprintf(atf, "%f %f %f %f\n", kweights[0], kweights[1], kweights[2], kweights[3]);
  
  int nr_local_ks = dftbp_get_nr_local_ks(&calculator);
  fprintf(atf, "nr_local_ks       :integer:0:\n%d\n", nr_local_ks); 

  // TODO use  typedef
  dftbp_register_dm_callback(&calculator, (void*)dm_callback, 0);
  dftbp_register_s_callback(&calculator, (void*)s_callback, 0);
  dftbp_register_h_callback(&calculator, (void*)h_callback, 0);
  
  assert(BASIS_SIZE == basis_size);
  assert(N_KPTS == 4);
  assert(N_SPIN == 1);

  /* Evaluate energy */
  double mermin_energy;
  dftbp_get_energy(&calculator, &mermin_energy);
  fprintf(atf, "mermin_energy       :real:0:\n%f\n",mermin_energy);
  
  
  for (int iK = 1; iK <= N_KPTS; ++iK)
  for (int iS = 1; iS <= N_SPIN; ++iS)
  {
    fprintf(atf, "dm_%d_%d       :complex:2:%d,%d\n",iK, iS,BASIS_SIZE, BASIS_SIZE);
    print_matrix(atf, &(dm[iK-1][iS-1][0][0]));
  }

  fprintf(atf, "overlap       :complex:2:%d,%d\n",BASIS_SIZE, BASIS_SIZE);
  print_matrix(atf, &(overlap[0][0]));

  fprintf(atf, "hamiltonian       :complex:2:%d,%d\n",BASIS_SIZE, BASIS_SIZE);
  print_matrix(atf, &(hamiltonian[0][0]));

  fclose(atf);
  return 0;
}
