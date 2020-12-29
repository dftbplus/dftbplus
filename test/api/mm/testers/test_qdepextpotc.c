/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2020  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

/**
 * Demonstrates the API with population dependant external potentials.
 *
 * Use it with the input in the qdepextpot/ folder.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dftbplus.h"
#include "testhelpers.h"

#define NR_QM_ATOMS 3
#define NR_MM_ATOMS 2

/**
 * Dummy helper for distance calculation
 */
double distance(const double *aa, const double *bb)
{
  return sqrt((aa[0] - bb[0]) * (aa[0] - bb[0])
              + (aa[1] - bb[1]) * (aa[1] - bb[1])
              + (aa[2] - bb[2]) * (aa[2] - bb[2]));
}

/**
 * Dummy helper to calculate an external potential caused by MM-charges.
 */
void calc_external_potential(int qmatoms, const double *qmcoords, int mmatoms,
                             const double *mmcoords, const double *mmcharges, double *extpot)
{
  const double *qmpos, *mmpos;
  double mmcharge, dist;
  int iatqm, iatmm;

  for (iatqm = 0; iatqm < qmatoms; ++iatqm) {
    extpot[iatqm] = 0.0;
    qmpos = &(qmcoords[3 * iatqm]);
    for (iatmm = 0; iatmm < mmatoms; ++iatmm) {
      mmpos = &(mmcoords[3 * iatmm]);
      mmcharge = mmcharges[iatmm];
      dist = distance(qmpos, mmpos);
      extpot[iatqm] += -mmcharge / dist;
    }
  }
}


/**
   Dummy helper routine to calculate the gradient of an external potential caused by MM-charges.
*/
void calc_external_potential_grad(int qmatoms, const double *qmcoords, int mmatoms,
                                  const double *mmcoords, const double *mmcharges,
                                  double *extpotgrad)
{
  const double *qmpos, *mmpos;
  double mmcharge, dist, dist3;
  int iatqm, iatmm;

  for (iatqm = 0; iatqm < qmatoms; ++iatqm) {
    qmpos = &(qmcoords[3 * iatqm]);
    extpotgrad[3 * iatqm + 0] = 0.0;
    extpotgrad[3 * iatqm + 1] = 0.0;
    extpotgrad[3 * iatqm + 2] = 0.0;
    for (iatmm = 0; iatmm < mmatoms; ++iatmm) {
      mmpos = &(mmcoords[3 * iatmm]);
      mmcharge = mmcharges[iatmm];
      dist = distance(qmpos, mmpos);
      dist3 = dist * dist * dist;
      extpotgrad[3 * iatqm + 0] += -mmcharge * (mmpos[0] - qmpos[0]) / dist3;
      extpotgrad[3 * iatqm + 1] += -mmcharge * (mmpos[1] - qmpos[1]) / dist3;
      extpotgrad[3 * iatqm + 2] += -mmcharge * (mmpos[2] - qmpos[2]) / dist3;
    }
  }
}


/**
 * Structure containing all the simulation data
 */
typedef struct {
  int mmatoms;
  int qmatoms;
  double *mmcoords;
  double *mmcharges;
  double *qmcoords;
  double *extpot;
  double *extpotgrad;
} Context;


/**
 * Call-back function calculating the external field caused by MM-charges.
 *
 * Note: This function has the signature of a generator for population-dependant external potential
 * only in order to demonstrate and test the usage of such callback functions. The function itself
 * does not make use of the population, and it that case it would be more efficient to pass the
 * population independant external potential via the dftbp_set_external_pot() function. Use this
 * kind of callback registration only, if you have an external potential (e.g. Kloppmann-Ohno),
 * which needs the populations on the atoms, and must be, therefore, updated in every SCC-cycle.
 *
 * Note: DFTB+ expects an external potential for the electrons. If you calculate an electrostatic
 * potential (electron having a negative charge), you must invert its sign.
 */
void get_external_potential(void *refptr, double *dqatom, double *extpotatom)
{
  Context *cont;

  cont = (Context *) refptr;
  calc_external_potential(cont->qmatoms, cont->qmcoords, cont->mmatoms, cont->mmcoords,
                          cont->mmcharges, extpotatom);
}


/**
 * Call-back function calculating the gradient of the external field caused by MM-charges.
 *
 * Note: This function has the signature of a generator for population-dependant external potential
 * only in order to demonstrate and test the usage of such callback functions. The function itself
 * does not make use of the population, and it that case it would be more efficient to pass the
 * population independant external potential via the dftbp_set_external_pot() function. Use this
 * kind of callback registration only, if you have an external potential (e.g. Kloppmann-Ohno),
 * which needs the populations on the atoms, and must be, therefore, updated in every SCC-cycle.
 *
 * Note: DFTB+ expects an external potential for the electrons. If you calculate an electrostatic
 * potential (electron having a negative charge), you must invert its sign.
 */
void get_external_potential_grad(void *refptr, double *dqatom, double *extpotgrad)
{
  Context *cont;

  cont = (Context *) refptr;
  calc_external_potential_grad(cont->qmatoms, cont->qmcoords, cont->mmatoms, cont->mmcoords,
                               cont->mmcharges, extpotgrad);
}


/**
 * Constructs the calculation context.
 *
 * This will be a simple type containing all the relevant variables which the callback function
 * has to deal with in order to calculate the external potential or its gradient.
 * Everything in atomic units.
 */
void initialize_context(Context *cont)
{
  const double MM_COORDS[3 * NR_MM_ATOMS] = {
    -0.944863438887178, -9.44863438887178, 1.70075418999692,
     4.34637181888102,  -5.85815332110050, 2.64561762888410
  };

  const double MM_CHARGES[2] = {2.5, -1.9};

  const double QM_COORDS[3 * NR_QM_ATOMS] = {
    0.0000000000000000, -1.8897259885789233, 0.0000000000000000,
    0.0000000000000000,  0.0000000000000000, 1.4797763915205659,
    0.0000000000000000,  0.0000000000000000, -1.4797763915205659
  };

  /* The potential of the first external charge will be calculated via the population-dependant
     potential generator (for demonstration purposes), while the second one will be stored and set
     as a simple static external potential.
  */
  cont->mmatoms = 1;
  cont->mmcoords = (double *) malloc(3 * cont->mmatoms * sizeof(double));
  memcpy(cont->mmcoords, MM_COORDS, 3 * cont->mmatoms * sizeof(double));
  cont->mmcharges = (double *) malloc(cont->mmatoms * sizeof(double));
  memcpy(cont->mmcharges, MM_CHARGES, cont->mmatoms * sizeof(double));

  cont->extpot = (double *) malloc(NR_QM_ATOMS * sizeof(double));
  calc_external_potential(
      NR_QM_ATOMS, QM_COORDS, (NR_MM_ATOMS - cont->mmatoms), MM_COORDS + 3 * cont->mmatoms,
      MM_CHARGES + cont->mmatoms, cont->extpot);

  cont->extpotgrad = (double *) malloc(3 * NR_QM_ATOMS * sizeof(double));
  calc_external_potential_grad(
      NR_QM_ATOMS, QM_COORDS, (NR_MM_ATOMS - cont->mmatoms), MM_COORDS + 3 * cont->mmatoms,
      MM_CHARGES + cont->mmatoms, cont->extpotgrad);

  cont->qmatoms = NR_QM_ATOMS;
  cont->qmcoords = (double *) malloc(3 * cont->qmatoms * sizeof(double));
  memcpy(cont->qmcoords, QM_COORDS, 3 * cont->qmatoms * sizeof(double));
}


/**
 * Destroys the calculational context.
 */
void finalize_context(Context *cont)
{
  free(cont->qmcoords);
  free(cont->mmcharges);
  free(cont->mmcoords);
  free(cont->extpot);
  free(cont->extpotgrad);
}


int main()
{
  DftbPlus calculator;
  DftbPlusInput input;
  DftbPlusAtomList dummyAtomList;

  Context cont;
  double mermin_energy;
  double *gradients, *charges;

  dummyAtomList.pDftbPlusAtomList = NULL;

  /* Fill up the context with all the relevant data */
  initialize_context(&cont);

  /* Initialize the calculator */
  dftbp_init(&calculator, NULL);
  /*dftbp_init(&calculator, "/dev/null");*/

  /* Parse the input file and store the input-tree */
  dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);

  /* Set up the calculator by processing the input tree */
  dftbp_process_input(&calculator, &input, &dummyAtomList);

  /* Register the callback functions calculating population dependent external potential */
  dftbp_register_ext_pot_generator(&calculator, &cont, get_external_potential,
                                   get_external_potential_grad);

  /* Send the current coordinates to deal with */
  dftbp_set_coords(&calculator, cont.qmcoords);

  /* Set current population-independent external field */
  dftbp_set_external_potential(&calculator, cont.extpot, cont.extpotgrad);

  /* Query the Mermin energy of the system */
  dftbp_get_energy(&calculator, &mermin_energy);
  printf("Expected Mermin-energy: %15.10f\n", -3.9854803392);
  printf("Obtained Mermin-energy: %15.10f\n", mermin_energy);

  /* Query the gradients on the atoms */
  gradients = (double *) malloc(3 * NR_QM_ATOMS * sizeof(double));
  dftbp_get_gradients(&calculator, gradients);
  printf("Expected gradient of atom 1: %15.10f %15.10f %15.10f\n", 0.017651363773,
         -0.183137129803, 0.003198251627);
  printf("Obtained gradient of atom 1: %15.10f %15.10f %15.10f\n", gradients[0], gradients[1],
         gradients[2]);

  /* Query the charges (electron population) on the atoms */
  charges = (double *) malloc(NR_QM_ATOMS * sizeof(double));
  dftbp_get_gross_charges(&calculator, charges);
  printf("Expected gross atom charges: %15.10f %15.10f %15.10f\n", -0.4943983279018598,
         0.264172212777779, 0.23022611512408198);
  printf("Obtained gross atom charges: %15.10f %15.10f %15.10f\n", charges[0], charges[1],
         charges[2]);


  /* Finalize the calculator */
  dftbp_final(&calculator);

  /* Destruct the calculational context */
  finalize_context(&cont);


  /* Save some data for the internal test system */
  dftbp_write_autotest_tag(NR_QM_ATOMS, 0, mermin_energy, gradients, NULL, charges, NULL);

  free(gradients);
  free(charges);

  return 0;
}
