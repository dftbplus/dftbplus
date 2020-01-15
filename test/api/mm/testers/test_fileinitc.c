/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2020  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "dftbplus.h"
#include "testhelpers.h"

/* This should match the number of atoms in the dftb_in.hsd */
#define NR_OF_ATOMS_SI2 2
#define NR_OF_ATOMS_H2O 3
#define MAX_ATOMS 3

#define NR_ITER 3

/*
  Reads in a dftb_in.hsd file then overrides some of its supplied setting, before calculating
  properties.
*/

void init_collective_variables(double *mermin_energy_total, double **gradients_total,
                               double **gross_charges_total)
{
  int ii;

  *gradients_total = (double *) calloc(3 * MAX_ATOMS, sizeof(double));
  *gross_charges_total = (double *) calloc(MAX_ATOMS, sizeof(double));
  *mermin_energy_total = 0.0;
  for (ii = 0; ii < MAX_ATOMS; ++ii) {
    (*gradients_total)[3 * ii] = 0.0;
    (*gradients_total)[3 * ii + 1] = 0.0;
    (*gradients_total)[3 * ii + 2] = 0.0;
    (*gross_charges_total)[ii] = 0.0;
  }
}


void update_collective_variables(int natom, double mermin_energy, double *gradients,
    double *gross_charges, double *mermin_energy_total, double *gradients_total,
    double *gross_charges_total)
{
  int ii;
  *mermin_energy_total += mermin_energy;
  for (ii = 0; ii < 3 * natom; ++ii) {
    gradients_total[ii] += gradients[ii];
  }
  for (ii = 0; ii < natom; ++ii) {
    gross_charges_total[ii] += gross_charges[ii];
  }
}


int main()
{
  DftbPlus calculator;
  DftbPlusInput input;

  /* Coordinates in row major format, atomic units */
  double coords_si2[3 * NR_OF_ATOMS_SI2] = {
    0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
    2.2639291987021915, 2.4639291987021915, 2.5639291987021915
  };


  /* Coordinates in row major format, atomic units */
  double coords_h2o[3 * NR_OF_ATOMS_H2O] = {
    0.00000000000E+00,  -0.10000000000E+01,   0.00000000000E+00,
    0.00000000000E+00,   0.00000000000E+00,   0.88306400000E+00,
    0.00000000000E+00,   0.00000000000E+00,  -0.78306400000E+00
  };

  /* Lattice vectors in row major format, atomic units */
  double latvecs_si2[3 * 3] = {
    5.2278583974043830, 5.1278583974043830, 0.0000000000000000,
    0.0000000000000000, 5.3278583974043830, 5.1278583974043830,
    5.1278583974043830, 0.0000000000000000, 5.4278583974043830
  };

  double mermin_energy, mermin_energy_total;
  double *gradients, *gradients_total, *gross_charges, *gross_charges_total;
  int natom, natom0, natom_total;
  int si2;

  /* Collective variables will hold the summed up results of multiple test runs */
  init_collective_variables(&mermin_energy_total, &gradients_total, &gross_charges_total);

  /* Dummy loop to test subsequent initialisations / finalisations of the DFTB+ object */
  for (int ii = 0; ii < NR_ITER; ++ii) {

    /* Use input for Si2 and H2O alternatingly, starting with Si2 */
    si2 = !(ii % 2);
    if (si2) {
      natom0 = NR_OF_ATOMS_SI2;
    } else {
      natom0 = NR_OF_ATOMS_H2O;
    }

    /* Initialise DFTB+ input tree from input in external file */
    dftbp_init(&calculator, NULL);
    if (si2) {
      dftbp_get_input_from_file(&calculator, "dftb_in.Si2.hsd", &input);
    } else {
      dftbp_get_input_from_file(&calculator, "dftb_in.H2O.hsd", &input);
    }
    dftbp_process_input(&calculator, &input);

    /* Check whether the calculator was initialized with the correct nr. of atoms */
    natom = dftbp_get_nr_atoms(&calculator);
    printf("Obtained nr. of atoms: %d\n", natom);
    printf("Expected nr. of atoms: %d\n", natom0);

    /* Override coordinates in the tree*/
    if (si2) {
      dftbp_set_coords_and_lattice_vecs(&calculator, coords_si2, latvecs_si2);
    } else {
      dftbp_set_coords(&calculator, coords_h2o);
    }

    /* evaluate energy */
    dftbp_get_energy(&calculator, &mermin_energy);
    printf("Obtained Mermin-energy: %15.10f\n", mermin_energy);
    if (si2) {
      printf("Expected Mermin-energy: %15.10f\n", -2.5933460731);
    } else {
      /* Energy not corrected yet! */
      printf("Expected Mermin-energy: %15.10f\n", -2.5933460731);
    }

    /* and gradients */
    gradients = (double *) calloc(3 * natom0, sizeof(double));
    dftbp_get_gradients(&calculator, gradients);
    printf("Obtained gradient of atom 1: %15.10f %15.10f %15.10f\n", gradients[3 * 0 + 0],
           gradients[3 * 0 + 1], gradients[3 * 0 + 2]);
    if (si2) {
      printf("Expected gradient of atom 1: %15.10f %15.10f %15.10f\n", -0.010321385989,
             -0.010321385989, -0.010321385989);
    } else { /* Needs to be corrected */
      printf("Expected gradient of atom 1: %15.10f %15.10f %15.10f\n", -0.010321385989,
             -0.010321385989, -0.010321385989);
    }

    /* and gross charges */
    gross_charges = (double *) calloc(natom0, sizeof(double));
    dftbp_get_gross_charges(&calculator, gross_charges);
    if (si2) {
      printf("Obtained gross charges: %15.10f %15.10f\n", gross_charges[0], gross_charges[1]);
    } else {
      printf("Obtained gross charges: %15.10f %15.10f %15.10f\n", gross_charges[0],
             gross_charges[1], gross_charges[2]);
    }

    update_collective_variables(natom0, mermin_energy, gradients, gross_charges,
                                &mermin_energy_total, gradients_total, gross_charges_total);

    /*  clean up */
    dftbp_final(&calculator);

    if (ii == NR_ITER - 1) {
      /* Save some data for the internal test system */
      dftbp_write_autotest_tag(MAX_ATOMS, 0, mermin_energy_total, gradients_total,
                               gross_charges_total, NULL);
    }

    free(gradients);
    free(gross_charges);
  }

  free(gradients_total);
  free(gross_charges_total);

  return 0;
}
