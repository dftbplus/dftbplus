/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2023  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

/**
 * MPI parallel calculation via the API and checks the result
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "dftbplus.h"
#include "testhelpers.h"

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("Hello world from processor %s, rank %d out of %d processors\n",
	 processor_name, world_rank, world_size);

  int major, minor, patch;
  dftbp_api(&major, &minor, &patch);
  if (world_rank == 0) {
    printf("DFTB+ API version %d.%d.%d\n", major, minor, patch);
  }

  // Set up an empty calculator
  DftbPlus calculator;
  dftbp_init_mpi(&calculator, NULL, MPI_Comm_c2f(MPI_COMM_WORLD));

  // Symmetrically initialise from an input file
  DftbPlusInput input;
  dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);
  dftbp_process_input(&calculator, &input);
  dftbp_input_final(&input);

  // Evaluate energy and forces
  double mermin_energy;
  dftbp_get_energy(&calculator, &mermin_energy);

  if (world_rank == 0) {
    printf("\nMermin free energy: %f\n", mermin_energy);
  }

  int nAtom = dftbp_get_nr_atoms(&calculator);

  if (world_rank == 0) {
    printf("\nAtoms: %i\n", nAtom);
  }

  /* Transfer 'nuclear' charge between some atoms at the start and end
     of the structure (i.e. some sort of virtual doping') */
  int nDelta = 3;
  if (nAtom > 2*nDelta) {

    // Change in reference valence charge in region
    double delta_z = 0.01;

    double* z0 = (double *) malloc(nAtom * sizeof(double));

    dftbp_get_ref_charges(&calculator, z0);

    if (world_rank == 0) {
      printf("Reference atomic occupations (q0)\n");
      for (int i=0; i<nAtom; i++) {
        printf("%i %f\n", i, z0[i]);
      }
    }

    for (int i=0; i<3; i++) {
      z0[i] += delta_z / (double) nDelta;
      z0[nAtom-i-1] -= delta_z / (double) nDelta;
    }

    if (world_rank == 0) {
      printf("Modified atomic occupations (q0)\n");
      for (int i=0; i<nAtom; i++) {
        printf("%i %f\n", i, z0[i]);
      }
    }

    // set charges in DFTB+
    dftbp_set_ref_charges(&calculator, z0);

    // Replacement energy in modified system
    dftbp_get_energy(&calculator, &mermin_energy);

    if (world_rank == 0) {
      printf("\nMermin free energy: %f\n", mermin_energy);
    }

    free(z0);
  }

  // Energy gradients (-forces0
  double* gradients = (double *) malloc(3 * nAtom * sizeof(double));
  dftbp_get_gradients(&calculator, gradients);

  // Gross charges of atoms
  double* charges = (double *) malloc(nAtom * sizeof(double));
  dftbp_get_gross_charges(&calculator, charges);

  // Clean up
  dftbp_final(&calculator);

  if (world_rank == 0) {
    /* Save some data for the internal test system */
    dftbp_write_autotest_tag(nAtom, 0, 0, mermin_energy, gradients, NULL, charges, NULL,
                             NULL, NULL);
  }

  free(gradients);
  free(charges);

  // Finalize the MPI environment.
  MPI_Finalize();
}
