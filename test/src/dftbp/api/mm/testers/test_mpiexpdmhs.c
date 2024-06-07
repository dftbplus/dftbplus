/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2024  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dftbplus.h"

#include "blacsutils.h"

int mpi_provided_threading, world_size, world_rank;

int basis_size = -1;
int n_spin = -1;
int n_kpts = -1;

typedef double hs_type;

hs_type *dm=0, *overlap=0, *hamiltonian=0;

void gather(int *blacs_descr, void *blacs_data, void *dest) {
  int blacs_ctx = blacs_descr[CTXT_];
  int sys_ctx = get_system_context(blacs_ctx);

  int gatherer_blacs_ctx = make_blacs_context(sys_ctx, 1, 1);
  int gathered_desc[DLEN_];
  blacs_desc_init(blacs_descr[M_], blacs_descr[N_], gatherer_blacs_ctx,
                  gathered_desc);

  int nprow, npcol, myrow, mycol;
  blacs_gridinfo_(&blacs_ctx, &nprow, &npcol, &myrow, &mycol);
  assert((gathered_desc[CTXT_] != -1) == ((myrow == 0) && (mycol == 0)));

  int ONE = 1;
  pdgemr2d_(&blacs_descr[M_], &blacs_descr[N_], blacs_data, &ONE, &ONE,
            blacs_descr, dest, &ONE, &ONE, gathered_desc, &blacs_ctx);
}

void dm_callback(void *aux_ptr, int iK, int iS, int *blacs_descr,
                 void *blacs_data, ASI_matrix_descr_t *matrix_descr) {
  if (matrix_descr->storage_type!=ASI_STORAGE_TYPE_DENSE_FULL)
  {
    abort();
  }
  memcpy(dm, blacs_data, sizeof(hs_type) * basis_size * basis_size);
  gather(blacs_descr, blacs_data, dm);
}

void s_callback(void *aux_ptr, int iK, int iS, int *blacs_descr,
                void *blacs_data, ASI_matrix_descr_t *matrix_descr) {
  if (matrix_descr->storage_type!=ASI_STORAGE_TYPE_DENSE_FULL)
  {
    abort();
  }
  memcpy(overlap, blacs_data, sizeof(hs_type) * basis_size * basis_size);
  gather(blacs_descr, blacs_data, overlap);
}

void h_callback(void *aux_ptr, int iK, int iS, int *blacs_descr,
                void *blacs_data, ASI_matrix_descr_t *matrix_descr) {
  if (matrix_descr->storage_type!=ASI_STORAGE_TYPE_DENSE_FULL)
  {
    abort();
  }
  memcpy(hamiltonian, blacs_data, sizeof(hs_type) * basis_size * basis_size);
  gather(blacs_descr, blacs_data, hamiltonian);
}

int main(int argc, char *argv[]) {
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,
                  &mpi_provided_threading);  // instead of MPI_Init
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  FILE *atf;
  if (world_rank == 0) {
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
  dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);
  dftbp_process_input(&calculator, &input);
  blacs_pinfo_(&mypnum, &nprocs);

  int natom = dftbp_get_nr_atoms(&calculator);
  basis_size = dftbp_get_basis_size(&calculator);
  _Bool is_hs_real = dftbp_is_hs_real(&calculator);
  n_spin = dftbp_get_nr_spin(&calculator);
  n_kpts = dftbp_nr_kpoints(&calculator);
  int nr_local_ks = dftbp_get_nr_local_ks(&calculator);
  int local_ks[n_kpts * n_spin];
  dftbp_get_local_ks(&calculator, local_ks);

  assert(is_hs_real);

  dm = (hs_type*) malloc(basis_size * basis_size * sizeof(hs_type));
  overlap = (hs_type*) malloc(basis_size * basis_size * sizeof(hs_type));
  hamiltonian = (hs_type*) malloc(basis_size * basis_size * sizeof(hs_type));


  dftbp_register_dm_callback(&calculator, dm_callback, 0);
  dftbp_register_s_callback(&calculator, s_callback, 0);
  dftbp_register_h_callback(&calculator, h_callback, 0);

  /* Evaluate energy */
  double mermin_energy;
  dftbp_get_energy(&calculator, &mermin_energy);

  double n_el = 0;
  for (int i = 0; i < basis_size; ++i)
    for (int j = 0; j < basis_size; ++j) {
      hs_type d = dm[j*basis_size+i];
      n_el += overlap[i*basis_size+j] * d;
    }

  if (world_rank == 0) {
    fprintf(atf, "natom          :integer:0:\n%d\n", natom);
    fprintf(atf, "basis_size     :integer:0:\n%d\n", basis_size);
    fprintf(atf, "is_hs_real     :integer:0:\n%d\n", (is_hs_real ? 1 : 0));
    fprintf(atf, "n_spin         :integer:0:\n%d\n", n_spin);
    fprintf(atf, "n_kpts         :integer:0:\n%d\n", n_kpts);
    fprintf(atf, "nr_local_ks    :integer:0:\n%d\n", nr_local_ks);
    fprintf(atf, "mermin_energy  :real:0:\n%f\n", mermin_energy);
    fprintf(atf, "n_el           :real:0:\n%f\n", n_el);

    fprintf(atf, "local_ks       :integer:1:%d\n", nr_local_ks * 2);
    for (size_t i = 0; i < nr_local_ks; ++i) {
      fprintf(atf, "%d %d", local_ks[i * 2], local_ks[i * 2 + 1]);
    }
    fprintf(atf, "\n");
    fclose(atf);
  }

  MPI_Finalize();
  return 0;
}
