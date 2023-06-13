/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2023  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dftbplus.h"
#include "testhelpers.h"

typedef double complex hsdC_type;
typedef double hsdR_type;

/* Need to make this all file scope global due to being a stateless
   protocol on this side. Arrays will be allocated to hold
   [n_kpts][n_spin][n_basis][n_basis] of data*/
static hsdC_type* dmC; // density matrix
static hsdC_type* overlapC;
static hsdC_type* hamiltonianC;
static hsdR_type* dmR;
static hsdR_type* overlapR;
static hsdR_type* hamiltonianR;
// Matrix properties
static _Bool is_real;

// As we already have global variables, store sizings here as well
static int n_kpts;
static int n_spin;
static int n_basis;

// Allocate contiguous memory for a complex array to hold
// [n_kpts][n_spin][n_basis][n_basis]
static void allocateCArray(hsdC_type** arrayRef) {
  int totalElements =  n_kpts * n_spin * n_basis * n_basis;
  *arrayRef = (hsdC_type*)malloc(totalElements * sizeof(hsdC_type));
  assert(*arrayRef != NULL);
}
// Allocate contiguous memory for a real array
static void allocateRArray(hsdR_type** arrayRef) {
  int totalElements =  n_kpts * n_spin * n_basis * n_basis;
  *arrayRef = (hsdR_type*)malloc(totalElements * sizeof(hsdR_type));
  assert(*arrayRef != NULL);
}

// Function to locate the elements location of the 4D array using
// [n_kpts][n_spin][n_basis][n_basis] assumed shape and Fortran
// indexing from 1 on the first two variables
inline static int findElement(int iKpt, int iSpin, int iBasis, int jBasis) {
  int i = iKpt -1;
  int j = iSpin -1;
  return (((i*n_spin) + j)*n_basis + iBasis)*n_basis + jBasis;
}

// Callback for density matrix access under ASI protocol, note iS and
// iK use Fortran array fencepost convention
void dm_callback(void *aux_ptr, int iK, int iS, int *blacs_descr, void *blacs_data) {
  if (is_real) {
    memcpy(&dmR[findElement(iK,iS,0,0)], blacs_data, sizeof(hsdR_type) * n_basis * n_basis);
  } else {
    memcpy(&dmC[findElement(iK,iS,0,0)], blacs_data, sizeof(hsdC_type) * n_basis * n_basis);
  }
}

// Callback for overlap matrix access
void s_callback(void *aux_ptr, int iK, int iS, int *blacs_descr, void *blacs_data) {
  if (is_real) {
    memcpy(&overlapR[findElement(iK,iS,0,0)], blacs_data, sizeof(hsdR_type) * n_basis * n_basis);
  } else {
    memcpy(&overlapC[findElement(iK,iS,0,0)], blacs_data, sizeof(hsdC_type) * n_basis * n_basis);
  }
}

// Callback for hamiltonian matrix access
void h_callback(void *aux_ptr, int iK, int iS, int *blacs_descr, void *blacs_data) {
  if (is_real) {
    memcpy(&hamiltonianR[findElement(iK,iS,0,0)], blacs_data, sizeof(hsdR_type) * n_basis * n_basis);
  } else {
    memcpy(&hamiltonianC[findElement(iK,iS,0,0)], blacs_data, sizeof(hsdC_type) * n_basis * n_basis);
  }
}

// Complex matrix print
static void print_Cmatrix(FILE *f, const hsdC_type *m) {
  for (int i = 0; i < n_basis; ++i) {
    // Upper triangle copied from lower
    for (int j = 0; j < i; ++j) {
      const hsdC_type dc = conj(m[j * n_basis + i]);
      fprintf(f, "%20.12e %20.12e\n", creal(dc), cimag(dc));
    }
    // Lower triangle
    for (int j = i; j < n_basis; ++j) {
      const hsdC_type dc = m[i * n_basis + j];
      fprintf(f, "%20.12e %20.12e\n", creal(dc), cimag(dc));
    }
  }
}

// Real matrix print
static void print_Rmatrix(FILE *f, const hsdR_type *m) {
  for (int i = 0; i < n_basis; ++i) {
    for (int j = 0; j < i; ++j) {
      // Upper triangle copied from lower
      fprintf(f, "%20.12e\n", m[j * n_basis + i]);
    }
    for (int j = i; j < n_basis; ++j) {
      fprintf(f, "%20.12e\n", m[i * n_basis + j]);
    }
  }
}


// Print out interface information
static void dftb_api_info() {

  int major, minor, patch;
  dftbp_api(&major, &minor, &patch);
  printf("API version %d.%d.%d\n", major, minor, patch);

  _Bool instsafe = dftbp_is_instance_safe();
  printf(instsafe ? "API is instance safe\n" : "API is NOT instance safe\n");

}

// clean up large stored matrices
static void cleanup() {
  if (is_real) {
    free(dmR);
    free(overlapR);
    free(hamiltonianR);
  } else {
    free(dmC);
    free(overlapC);
    free(hamiltonianC);
  }
}

int main() {
  FILE *atf = fopen("autotest.tag", "w+");

  DftbPlus calculator;
  DftbPlusInput input;

  dftb_api_info();

  /* Initialize DFTB+ input tree from input in external file */
  dftbp_init(&calculator, NULL);
  dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);
  dftbp_process_input(&calculator, &input);

  int natom = dftbp_get_nr_atoms(&calculator);
  fprintf(atf, "natom       :integer:0:\n%d\n", natom);

  is_real = dftbp_is_hs_real(&calculator);
  fprintf(atf, "is_hs_real       :integer:0:\n%d\n", (is_real ? 1 : 0));

  n_basis = dftbp_get_basis_size(&calculator);
  fprintf(atf, "basis_size       :integer:0:\n%d\n", n_basis);

  n_spin = dftbp_get_nr_spin(&calculator);
  fprintf(atf, "n_spin       :integer:0:\n%d\n", n_spin);

  n_kpts = dftbp_nr_kpoints(&calculator);
  fprintf(atf, "n_kpts       :integer:0:\n%d\n", n_kpts);

  hsdR_type kweights[n_kpts];
  dftbp_get_kweights(&calculator, kweights);
  fprintf(atf, "kweights       :real:1:%d\n", n_kpts);
  for (int i = 0; i < n_kpts -1; ++i) {
    fprintf(atf, "%14.10f ", kweights[i]);
  }
  fprintf(atf, "%14.10f\n", kweights[n_kpts -1]);

  int nr_local_ks = dftbp_get_nr_local_ks(&calculator);
  fprintf(atf, "nr_local_ks       :integer:0:\n%d\n", nr_local_ks);

  int local_ks[nr_local_ks * 2];
  dftbp_get_local_ks(&calculator, local_ks);
  fprintf(atf, "local_ks       :integer:2:%d,%d\n", nr_local_ks, 2);
  for (int i = 0; i < nr_local_ks; ++i) {
    fprintf(atf, "%i %i\n", local_ks[i * 2], local_ks[i * 2 + 1]);
  }

  if (is_real) {
    allocateRArray(&dmR);
    allocateRArray(&overlapR);
    allocateRArray(&hamiltonianR);
  } else {
    allocateCArray(&dmC);
    allocateCArray(&overlapC);
    allocateCArray(&hamiltonianC);
  }

  dftbp_register_dm_callback(&calculator, dm_callback, 0);
  dftbp_register_s_callback(&calculator, s_callback, 0);
  dftbp_register_h_callback(&calculator, h_callback, 0);

  /* Evaluate energy */
  hsdR_type mermin_energy;
  dftbp_get_energy(&calculator, &mermin_energy);
  fprintf(atf, "mermin_energy       :real:0:\n%14.10f\n", mermin_energy);

  if (is_real) {
    hsdR_type Hsum = 0, Ssum = 0;

    for (int iK = 1; iK <= n_kpts; ++iK)
      for (int iS = 1; iS <= n_spin; ++iS) {
        fprintf(atf, "dm_%d_%d       :real:2:%d,%d\n", iK, iS, n_basis, n_basis);
        print_Rmatrix(atf, &dmR[findElement(iK, iS, 0, 0)]);

        fprintf(atf, "overlap_%d_%d       :real:2:%d,%d\n", iK, iS, n_basis, n_basis);
        print_Rmatrix(atf, &overlapR[findElement(iK, iS, 0, 0)]);

        fprintf(atf, "hamiltonian_%d_%d       :real:2:%d,%d\n", iK, iS, n_basis, n_basis);
        print_Rmatrix(atf, &hamiltonianR[findElement(iK, iS, 0, 0)]);

        for (int i = 0; i < n_basis; ++i) {
          hsdR_type d = dmR[findElement(iK, iS, i, i)];
          hsdR_type s = overlapR[findElement(iK, iS, i, i)];
          hsdR_type h = hamiltonianR[findElement(iK, iS, i, i)];
          Ssum += kweights[iK-1] * d * s;
          Hsum += kweights[iK-1] * d * h;
          for (int j = i+1; j < n_basis; ++j) {
            hsdR_type d = dmR[findElement(iK, iS, i, j)];
            hsdR_type s = overlapR[findElement(iK, iS, i, j)];
            hsdR_type h = hamiltonianR[findElement(iK, iS, i, j)];
            Ssum += kweights[iK-1] * 2.0 * d * s;
            Hsum += kweights[iK-1] * 2.0 * d * h;
          }
        }
      }
    fprintf(atf, "Ssum         :real:0:\n%14.10f\n", Ssum);
    fprintf(atf, "Hsum         :real:0:\n%14.10f\n", Hsum);

  } else {

    hsdC_type Hsum = 0, Ssum = 0;
    for (int iK = 1; iK <= n_kpts; ++iK)
      for (int iS = 1; iS <= n_spin; ++iS) {
        fprintf(atf, "dm_%d_%d       :complex:2:%d,%d\n", iK, iS, n_basis, n_basis);
        print_Cmatrix(atf, &dmC[findElement(iK, iS, 0, 0)]);

        fprintf(atf, "overlap_%d_%d       :complex:2:%d,%d\n", iK, iS, n_basis, n_basis);
        print_Cmatrix(atf, &overlapC[findElement(iK, iS, 0, 0)]);

        fprintf(atf, "hamiltonian_%d_%d       :complex:2:%d,%d\n", iK, iS, n_basis, n_basis);
        print_Cmatrix(atf, &hamiltonianC[findElement(iK, iS, 0, 0)]);

        for (int i = 0; i < n_basis; ++i) {
          hsdC_type d = dmC[findElement(iK, iS, i, i)];
          hsdC_type s = overlapC[findElement(iK, iS, i, i)];
          hsdC_type h = hamiltonianC[findElement(iK, iS, i, i)];
          Ssum += kweights[iK-1] * d * s;
          Hsum += kweights[iK-1] * d * h;
          for (int j = i+1; j < n_basis; ++j) {
            hsdC_type d = dmC[findElement(iK, iS, i, j)];
            hsdC_type s = overlapC[findElement(iK, iS, i, j)];
            hsdC_type h = hamiltonianC[findElement(iK, iS, i, j)];
            Ssum += kweights[iK-1] * (d * conj(s) + conj(d) * s);
            Hsum += kweights[iK-1] * (d * conj(h) + conj(d) * h);
          }
        }
      }
    fprintf(atf, "Ssum         :complex:0:\n%14.10f %14.10f\n", creal(Ssum), cimag(Ssum));
    fprintf(atf, "Hsum         :complex:0:\n%14.10f %14.10f\n", creal(Hsum), cimag(Hsum));
  }

  cleanup();

  fclose(atf);
  return 0;
}
