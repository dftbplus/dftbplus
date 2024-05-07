/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2024  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <dlfcn.h>

void *init_plugin(const char *filename) {
  void *handle = dlopen(filename, RTLD_NOW);
  if (handle == NULL) return NULL;

  void (*initfunc)();
  initfunc = dlsym(handle, "init");
  if (initfunc != NULL) {
    (*initfunc)();
  }

  return handle;
}

void final_plugin(void *handle) {
  void (*finalfunc)();
  finalfunc = dlsym(handle, "final");
  if (finalfunc != NULL) {
    (*finalfunc)();
  }

  dlclose(handle);
}

int provides_plugin(void *handle, const char *func) {
  void *funcpointer = dlsym(handle, func);
  return funcpointer == NULL ? 0 : 1;
}

int call_getSKIntegrals(void *handle, int nSk, double *sk, double dist, int atom1, int atom2,
    int sp1, int sp2) {
  void (*func)(int, double *, double, int, int, int, int);
  func = dlsym(handle, "getSKIntegrals");
  if (func != NULL) {
    (*func)(nSk, sk, dist, atom1, atom2, sp1, sp2);
    return 1;
  }

  return 0;
}

int call_setNeighbourList(void *handle, int nAtoms, double **coords, int *img2CentCell) {
  void (*func)(int, double *, int *);
  func = dlsym(handle, "setNeighbourList");
  if (func != NULL) {
    (*func)(nAtoms, coords, img2CentCell);
    return 1;
  }

  return 0;
}
