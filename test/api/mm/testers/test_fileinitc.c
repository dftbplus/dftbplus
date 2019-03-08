/**************************************************************************************************/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2018  DFTB+ developers group                                                    */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/**************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "dftbplus.h"
#include "testhelpers.h"

/* This should match the number of atoms in the dftb_in.hsd */
#define NR_OF_ATOMS 2

/*
Reads in a dftb_in.hsd file then overrides some of its supplied setting, before calculating
properties. You should provide the dftb_in.hsd and skfiles as found in the
test/prog/dftb+/non-scc/Si_2/ folder
*/
int main()
{
  DftbPlus calculator;
  DftbPlusInput input;

  /* Coordinates in row major format, atomic units */
  double coords[3 * NR_OF_ATOMS] = {
    0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
    2.5639291987021915, 2.5639291987021915, 2.5639291987021915
  };

  /* Lattice vectors in row major format, atomic units */
  double latvecs[3 * 3] = {
    5.1278583974043830, 5.1278583974043830, 0.0000000000000000,
    0.0000000000000000, 5.1278583974043830, 5.1278583974043830,
    5.1278583974043830, 0.0000000000000000, 5.1278583974043830
  };

  double mermin_energy;
  double *gradients, *gross_charges;
  int natom;

  /* Initialise DFTB+ input tree from input in external file */
  dftbp_init(&calculator, NULL);
  dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);
  dftbp_process_input(&calculator, &input);

  /* Check whether the calculator was initialized with the correct nr. of atoms */
  natom = dftbp_get_nr_atoms(&calculator);
  printf("Expected nr. of atoms: %d\n", NR_OF_ATOMS);
  printf("Obtained nr. of atoms: %d\n", natom);

  /* Override coordinates in the tree*/
  dftbp_set_coords_and_lattice_vecs(&calculator, coords, latvecs);

  /* evaluate energy */
  dftbp_get_energy(&calculator, &mermin_energy);
  printf("Expected Mermin-energy: %15.10f\n", -2.5933460731);

  /* and gradients */
  gradients = (double *) calloc(3 * NR_OF_ATOMS, sizeof(double));
  dftbp_get_gradients(&calculator, gradients);
  printf("Expected gradient of atom 1: %15.10f %15.10f %15.10f\n", -0.010321385989, -0.010321385989,
         -0.010321385989);

  /* and gross charges */
  gross_charges = (double *) calloc(NR_OF_ATOMS, sizeof(double));
  dftbp_get_gross_charges(&calculator, gross_charges);
  printf("Obtained gross charges: %15.10f %15.10f\n", gross_charges[0], gross_charges[1]);

  /*  clean up */
  dftbp_final(&calculator);


  /* Save some data for the internal test system */
  dftbp_write_autotest_tag(NR_OF_ATOMS, 0, mermin_energy, gradients, gross_charges, NULL);

  return 0;
}
