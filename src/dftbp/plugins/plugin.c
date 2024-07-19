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

int version_plugin(void *handle, int major, int minor) {
  int (*versionfunc)(int, int);
  versionfunc = dlsym(handle, "version");
  if (versionfunc != NULL) {
    return (*versionfunc)(major, minor);
  }

  return 0;
}

int provides_plugin(void *handle, const char *func) {
  void *funcpointer = dlsym(handle, func);
  return funcpointer == NULL ? 0 : 1;
}

int call_getSKIntegrals(void *handle, int nSkgrid, int nSkIntg, double *skTab, double dist,
    int atom1, int atom2, int species1, int species2, int HorS, double interdist) {
  int (*func)(int, int, double *, double, int, int, int, int, int, double);
  func = dlsym(handle, "getSKIntegrals");
  if (func != NULL) {
    return (*func)(nSkgrid, nSkIntg, skTab, dist, atom1, atom2, species1, species2, HorS, interdist);
  }

  return 0;
}

void call_setNeighbourList(void *handle, int nAtoms, int nAtomsCent, double *coords,
    int *img2CentCell, int *iNeighbour, double *neightDist2) {
  int (*func)(int, int, double *, int *, int *, double *);
  func = dlsym(handle, "setNeighbourList");
  if (func != NULL) {
    (*func)(nAtoms, nAtomsCent, coords, img2CentCell, iNeighbour, neightDist2);
  }
}
