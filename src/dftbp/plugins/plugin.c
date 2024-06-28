/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2024  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <dlfcn.h>

void *init_plugin(const char *filename) {
  void *handle = dlopen(filename, RTLD_NOW | RTLD_GLOBAL);
  if (handle == NULL) return NULL;

  int (*initfunc)();
  initfunc = dlsym(handle, "init");
  if (initfunc != NULL) {
    int success = (*initfunc)();
    if (success != 1) {
      return NULL;
    }
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
    int species1, int species2) {
  int (*func)(int, double *, double, int, int, int, int);
  func = dlsym(handle, "getSKIntegrals");
  if (func != NULL) {
    return (*func)(nSk, sk, dist, atom1, atom2, species1, species2);
  }

  return 0;
}

void call_setNeighbourList(void *handle, int nAtoms, double *coords, int *img2CentCell) {
  int (*func)(int, double *, int *);
  func = dlsym(handle, "setNeighbourList");
  if (func != NULL) {
    (*func)(nAtoms, coords, img2CentCell);
  }
}
