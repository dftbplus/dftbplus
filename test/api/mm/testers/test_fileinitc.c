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

/* This should match the number of atoms in dftb_in.hsd */
#define NR_OF_ATOMS_SI2 2
#define NR_OF_ATOMS_H2O 3
#define MAX_ATOMS 3

#define NR_ITER 3

/*
  Reads in a dftb_in.hsd file then overrides some of
  its supplied settings, before calculating properties.
*/

void init_collective_variables(double *mermin_energy_total, double **gradients_total,
                               double **stress_tensor_total, double **gross_charges_total)
{
  int ii;

  *gradients_total = (double *) calloc(3 * MAX_ATOMS, sizeof(double));
  *stress_tensor_total = (double *) calloc(9, sizeof(double));
  *gross_charges_total = (double *) calloc(MAX_ATOMS, sizeof(double));
  *mermin_energy_total = 0.0;
  for (ii = 0; ii < MAX_ATOMS; ++ii) {
    (*gradients_total)[3 * ii] = 0.0;
    (*gradients_total)[3 * ii + 1] = 0.0;
    (*gradients_total)[3 * ii + 2] = 0.0;
    (*gross_charges_total)[ii] = 0.0;
  }
  for (ii = 0; ii < 9; ++ii) {
    (*stress_tensor_total)[ii] = 0.0;
  }
}


void update_collective_variables(int natom, double mermin_energy, double *gradients,
    double *stress_tensor, double *gross_charges, double *mermin_energy_total, 
    double *gradients_total, double *stress_tensor_total, double *gross_charges_total)
{
  int ii;
  *mermin_energy_total += mermin_energy;
  for (ii = 0; ii < 3 * natom; ++ii) {
    gradients_total[ii] += gradients[ii];
  }
  for (ii = 0; ii < natom; ++ii) {
    gross_charges_total[ii] += gross_charges[ii];
  }
  for (ii = 0; ii < 9; ++ii)
  {
    stress_tensor_total[ii] += stress_tensor[ii];
  }
}


int main()
{
  DftbPlus calculator;
  DftbPlusInput input;
  DftbPlusAtomList dummyAtomList;

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

  /* Arbitrary coordinate origin in atomic units */
  double origin_si2[3] = {
    5.2278583974043830, 5.1278583974043830, 0.0000000000000000
  };

  double mermin_energy, mermin_energy_total;
  double *gradients, *gradients_total, *stress_tensor, *stress_tensor_total, *gross_charges, *gross_charges_total;
  int natom, natom0;
  int si2, ii, ij;
  int major, minor, patch;

  dummyAtomList.pDftbPlusAtomList = NULL;

  dftbp_api(&major, &minor, &patch);
  printf("API version %d.%d.%d\n", major, minor, patch);

  /* Collective variables will hold the summed up results of multiple test runs */
  init_collective_variables(&mermin_energy_total, &gradients_total, 
                         &stress_tensor_total, &gross_charges_total);

  /* Dummy loop to test subsequent initializations / finalizations of the DFTB+ object */
  for (ii = 0; ii < NR_ITER; ++ii) {

    /* Use input for Si2 and H2O alternatingly, starting with Si2 */
    si2 = !(ii % 2);
    if (si2) {
      natom0 = NR_OF_ATOMS_SI2;
    } else {
      natom0 = NR_OF_ATOMS_H2O;
    }

    /* Initialize DFTB+ input tree from input in external file */
    dftbp_init(&calculator, NULL);
    if (si2) {
      dftbp_get_input_from_file(&calculator, "dftb_in.Si2.hsd", &input);
    } else {
      dftbp_get_input_from_file(&calculator, "dftb_in.H2O.hsd", &input);
    }
    dftbp_process_input(&calculator, &input, &dummyAtomList);

    /* Check whether the calculator was initialized with the correct nr. of atoms */
    natom = dftbp_get_nr_atoms(&calculator);
    printf("Obtained nr. of atoms: %d\n", natom);
    printf("Expected nr. of atoms: %d\n", natom0);

    /* Override coordinates in the tree*/
    if (si2) {
      if (ii == 0) {
        dftbp_set_coords_and_lattice_vecs(&calculator, coords_si2, latvecs_si2);
      } else {
        /* Set an origin for a periodic cell - for most cases this is
           arbitrary and not needed */
        dftbp_set_coords_lattice_origin(&calculator, coords_si2, latvecs_si2, origin_si2);
      }

    } else {
      dftbp_set_coords(&calculator, coords_h2o);
    }

    /* Evaluate energy */
    dftbp_get_energy(&calculator, &mermin_energy);
    printf("Obtained Mermin-energy: %15.10f\n", mermin_energy);
    if (si2) {
      printf("Expected Mermin-energy: %15.10f\n", -2.5897497363);
    } else {
      printf("Expected Mermin-energy: %15.10f\n", -3.7534584715);
    }

    /* Evaluate gradients */
    gradients = (double *) calloc(3 * natom0, sizeof(double));
    dftbp_get_gradients(&calculator, gradients);
    printf("Obtained gradient of atom 1: %15.10f %15.10f %15.10f\n", gradients[3 * 0 + 0],
           gradients[3 * 0 + 1], gradients[3 * 0 + 2]);
    if (si2) {
      printf("Expected gradient of atom 1: %15.10f %15.10f %15.10f\n", 0.0306186399,
	     0.0026710677, -0.0007231241);
    } else {
      printf("Expected gradient of atom 1: %15.10f %15.10f %15.10f\n", 0.0178274222,
	     1.1052419875, -0.0835794936);
    }

    /* Evaluate stress tensor (if the system is periodic) */
    stress_tensor = (double *) calloc(9, sizeof(double));
    if(si2) {
      dftbp_get_stress_tensor(&calculator, stress_tensor);
      printf("Obtained diagonal elements of stress tensor: %15.10f %15.10f %15.10f\n", 
           stress_tensor[0], stress_tensor[4], stress_tensor[8]);
      printf("Expected diagonal elements of stress tensor: %15.10f %15.10f %15.10f\n", 
	     -0.0001163360, -0.0001080554, -0.0001614502);
    } else {
      for (ij = 0; ij < 9; ++ij) {
        stress_tensor[ij] = 0.0;
      }
    }

    /* Evaluate Gross charges */
    gross_charges = (double *) calloc(natom0, sizeof(double));
    dftbp_get_gross_charges(&calculator, gross_charges);
    if (si2) {
      printf("Obtained gross charges: %15.10f %15.10f\n", gross_charges[0], gross_charges[1]);
      printf("Expected gross charges: %15.10f %15.10f\n", 0.0000000000, 0.0000000000);
    } else {
      printf("Obtained gross charges: %15.10f %15.10f %15.10f\n", gross_charges[0],
             gross_charges[1], gross_charges[2]);
      printf("Expected gross charges: %15.10f %15.10f %15.10f\n", -0.6519945363,
	     0.3314953102, 0.3204992261);
    }

    update_collective_variables(natom0, mermin_energy, gradients, stress_tensor, 
                                gross_charges, &mermin_energy_total, gradients_total, 
                                stress_tensor_total, gross_charges_total);

    /*  Clean up */
    dftbp_final(&calculator);

    if (ii == NR_ITER - 1) {
      /*
	 Save some data in autotest.tag for the internal test system.
         In order to include the results of all performed calculations,
         the summed up values (collective variables) are written out.
      */
      dftbp_write_autotest_tag(MAX_ATOMS, 0, mermin_energy_total, gradients_total,
                               stress_tensor_total, gross_charges_total, NULL);
    }

    free(gradients);
    free(stress_tensor);
    free(gross_charges);
  }

  free(gradients_total);
  free(stress_tensor_total);
  free(gross_charges_total);

  return 0;
}
