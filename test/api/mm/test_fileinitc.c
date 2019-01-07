/**************************************************************************************************/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2018  DFTB+ developers group                                                    */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/**************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "dftbplus.h"

#define NR_OF_ATOMS 2


int main()
{
  DftbPlus calculator;
  DftbPlusInput input;

  /* Coordinates in row major format */
  double coords[3 * NR_OF_ATOMS] = {
    0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
    2.5639291987021915, 2.5639291987021915, 2.5639291987021915
  };

  /* Lattice vectors in row major format */
  double latvecs[3 * 3] = {
    5.1278583974043830, 5.1278583974043830, 0.0000000000000000,
    0.0000000000000000, 5.1278583974043830, 5.1278583974043830,
    5.1278583974043830, 0.0000000000000000, 5.1278583974043830
  };

  double mermin_energy;
  double *gradients;
  int natom;

  dftbp_init(&calculator, NULL);
  dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);
  dftbp_process_input(&calculator, &input);

  /* Check whether the calculator was initialized with the correct nr. of atoms */
  natom = dftbp_get_nr_atoms(&calculator);
  printf("Expected nr. of atoms: %d\n", NR_OF_ATOMS);
  printf("Obtained nr. of atoms: %d\n", natom);

  dftbp_set_coords_and_lattice_vecs(&calculator, coords, latvecs);

  dftbp_get_energy(&calculator, &mermin_energy);
  printf("Expected Mermin-energy: %15.10f\n", -2.5933460731);
  printf("Obtained Mermin-energy: %15.10f\n", mermin_energy);

  gradients = (double *) calloc(3 * NR_OF_ATOMS, sizeof(double));
  dftbp_get_gradients(&calculator, gradients);
  printf("Expected gradient of atom 1: %15.10f %15.10f %15.10f\n", -0.010321385989, -0.010321385989,
         -0.010321385989);
  printf("Obtained gradient of atom 1: %15.10f %15.10f %15.10f\n", gradients[0], gradients[1],
         gradients[2]);

  dftbp_final(&calculator);

  return 0;
}
