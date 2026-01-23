/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2025  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

/**
 * Demonstrates the API with derivative response calculations.
 *
 * Use it with the input in the test/src/dftbp/api/mm/testcases/deriv_responsesc/ folder.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dftbplus.h"
#include "testhelpers.h"

int main()
{
  DftbPlus calculator;
  DftbPlusInput input;

  double mermin_energy;
  bool isPeriodic;

  /* Initialize the calculator */
  dftbp_init(&calculator, NULL);
  /*dftbp_init(&calculator, "/dev/null");*/

  /* Parse the input file and store the input-tree */
  dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);

  /* Set up the calculator by processing the input tree */
  dftbp_process_input(&calculator, &input);
  dftbp_input_final(&input);

  int natom = dftbp_get_nr_atoms(&calculator);
  int nextchrg = dftbp_get_nr_extchrg(&calculator);

  double *dqdx;
  int nAtomsWRT;
  int *atomsWRT;
  double *dqdxExt;
  int nChargesWRT;
  int *chargesWRT;

  /* Query the Mermin energy of the system */
  dftbp_get_energy(&calculator, &mermin_energy);
  printf("Expected Mermin-energy: %15.10f\n", -3.9854803392);
  printf("Obtained Mermin-energy: %15.10f\n", mermin_energy);

  double *latvecs = (double *) malloc(3 * 3 * sizeof(double));
  dftbp_get_latvecs(&calculator, &isPeriodic, latvecs);

  printf("\nSystem with %i atoms, %i external charges\n\n", natom, nextchrg);

  if (isPeriodic) {

    printf("Lattice vectors (Bohr)\n");
    for (int ii=0; ii<3; ii++) {
      printf("%f %f %f\n", latvecs[3*ii], latvecs[3*ii+1], latvecs[3*ii+2]);
    }
    printf("\n");

    /* Periodic, so not calculating derivatives wrt internal atomic coordinates */
    nAtomsWRT = 0;
    atomsWRT = NULL;
    dqdx = NULL;

  } else {

    nAtomsWRT = natom;
    atomsWRT = (int *) malloc(nAtomsWRT * sizeof(int));
    for (int ii = 0; ii < nAtomsWRT; ii++) {
      atomsWRT[ii] = ii+1;
    }
    dqdx = (double *) malloc(natom * 3 * nAtomsWRT * sizeof(double));

  }

  if (nextchrg > 0) {

    nChargesWRT = nextchrg;
    chargesWRT = (int *) malloc(nChargesWRT * sizeof(int));
    for (int ii = 0; ii < nChargesWRT; ii++) {
      chargesWRT[ii] = ii+1;
    }
    dqdxExt = (double *) malloc(natom * 3 * nChargesWRT * sizeof(double));

  } else {

    /* No external charges, so no derivatives */
    nChargesWRT = 0;
    chargesWRT = NULL;
    dqdxExt = NULL;

  }

  printf("Evaluating response wrt %i atoms and %i external charges\n\n", nAtomsWRT, nChargesWRT);

  /* Response calculations */
  dftbp_get_charge_derivatives_select(&calculator, nAtomsWRT, atomsWRT, dqdx, nChargesWRT,
				      chargesWRT, dqdxExt);

  if (! isPeriodic) {

    printf("\n WRT atoms\n");
    for (int ii = 0; ii < nAtomsWRT; ii++) {
      printf("wrt atom %d\n", ii+1);
      for (int jj = 0; jj < 3; jj++) {
        for (int kk = 0; kk < natom; kk++) {
          printf(" %f", dqdx[kk+jj*natom+ii*(3*nAtomsWRT)]);
        }
        printf("\n");
      }
      printf("\n");
    }

  }

  if (nextchrg > 0) {

    printf("\n WRT external charges\n");
    for (int ii = 0; ii < nChargesWRT; ii++) {
      printf("wrt charge %d\n", ii+1);
      for (int jj = 0; jj < 3; jj++) {
        for (int kk = 0; kk < natom; kk++) {
          printf(" %f", dqdxExt[kk+jj*natom+ii*(3*natom)]);
        }
        printf("\n");
      }
      printf("\n");
    }

  }

  /* Finalize the calculator */
  dftbp_final(&calculator);

  /* ave some data for the internal test system */
  dftbp_write_autotest_tag(natom, nextchrg, 0, mermin_energy, NULL, NULL, NULL, NULL, NULL, NULL, dqdx, dqdxExt);

  return 0;
}
