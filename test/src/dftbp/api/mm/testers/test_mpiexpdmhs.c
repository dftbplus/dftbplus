/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2022  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>

#include "dftbplus.h"

int mpi_provided_threading, world_size, world_rank;

#define BASIS_SIZE 8
#define N_KPTS 1
#define N_SPIN 2

double dm[BASIS_SIZE][BASIS_SIZE];
double overlap[BASIS_SIZE][BASIS_SIZE];
double hamiltonian[BASIS_SIZE][BASIS_SIZE];

void dm_callback(void *aux_ptr, int iK, int iS, int *blacs_descr, void *blacs_data)
{
  memcpy(dm, blacs_data, sizeof(double)*BASIS_SIZE*BASIS_SIZE);
}

void s_callback(void *aux_ptr, int *blacs_descr, void *blacs_data)
{
  memcpy(overlap, blacs_data, sizeof(double)*BASIS_SIZE*BASIS_SIZE);
}

void h_callback(void *aux_ptr, int *blacs_descr, void *blacs_data)
{
  memcpy(hamiltonian, blacs_data, sizeof(double)*BASIS_SIZE*BASIS_SIZE);
}


int main(int argc, char *argv[])
{
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_threading); // instead of MPI_Init
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  
  FILE *atf;
  if (world_rank == 0)
  {
    atf = fopen("autotest.tag", "w+");
  }

  int mypnum, nprocs;
  blacs_pinfo_(&mypnum, &nprocs);
  printf("mypnum = %d nprocs = %d @ %d\n", mypnum, nprocs, world_rank);

  int major, minor, patch;
  dftbp_api(&major, &minor, &patch);
  printf("API version %d.%d.%d\n", major, minor, patch);

  _Bool instsafe = dftbp_is_instance_safe();
  printf(instsafe ? "API is instance safe\n" : "API is NOT instance safe\n");
  
  DftbPlus calculator;
  const MPI_Fint f_mpi_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
  printf("world_size=%d @ %d\n", world_size, world_rank);

  dftbp_init_mpi(&calculator, NULL, f_mpi_comm);
  DftbPlusInput input;
  dftbp_get_input_from_file(&calculator, "dftb_in.O2.hsd", &input);
  dftbp_process_input(&calculator, &input);
  blacs_pinfo_(&mypnum, &nprocs);

  int natom = dftbp_get_nr_atoms(&calculator);
  int basis_size = dftbp_get_basis_size(&calculator);
  _Bool is_hs_real = dftbp_is_hs_real(&calculator);
  int n_spin = dftbp_get_nr_spin(&calculator);
  int n_kpts = dftbp_nr_kpoints(&calculator);
  int nr_local_ks = dftbp_get_nr_local_ks(&calculator);
  int local_ks[N_KPTS * N_SPIN];
  dftbp_get_local_ks(&calculator, local_ks);
  
  assert(BASIS_SIZE == basis_size);
  assert(N_KPTS == n_kpts);
  assert(N_SPIN == n_spin);


  // TODO use  typedef
  dftbp_register_dm_callback(&calculator, (void*)dm_callback, 0);
  dftbp_register_s_callback(&calculator, (void*)s_callback, 0);
  dftbp_register_h_callback(&calculator, (void*)h_callback, 0);

  /* Evaluate energy */
  double mermin_energy;
  dftbp_get_energy(&calculator, &mermin_energy);

  double n_el = 0;
  for(int i = 0; i < BASIS_SIZE; ++i)
  for(int j = 0; j < BASIS_SIZE; ++j)
  {
    double d = dm[i][j];
    if (i != j)
      d += dm[j][i];
    n_el += overlap[i][j] * d;
  }

  if (world_rank == 0)
  {
    fprintf(atf, "natom          :integer:0:\n%d\n",natom);
    fprintf(atf, "basis_size     :integer:0:\n%d\n", basis_size);
    fprintf(atf, "is_hs_real     :integer:0:\n%d\n", (is_hs_real ? 1 : 0) );
    fprintf(atf, "n_spin         :integer:0:\n%d\n", n_spin);
    fprintf(atf, "n_kpts         :integer:0:\n%d\n", n_kpts);
    fprintf(atf, "nr_local_ks    :integer:0:\n%d\n", nr_local_ks); 
    fprintf(atf, "mermin_energy  :real:0:\n%f\n",mermin_energy);
    fprintf(atf, "n_el           :real:0:\n%f\n",n_el);
    
    fprintf(atf, "local_ks       :integer:1:%d\n", nr_local_ks*2);
    for (size_t i=0; i < nr_local_ks; ++i)
    {
      fprintf(atf, "%d %d", local_ks[i*2], local_ks[i*2 + 1]);
    }
    fprintf(atf, "\n");
    fclose(atf);
  }

  
  MPI_Finalize();
  return 0;
}
