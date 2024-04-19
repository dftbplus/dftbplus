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

int call_getSKIntegrals(void *handle, double *sk, double dist, int sp1, int sp2) {
  void (*func)(double *, double, int, int);
  func = dlsym(handle, "getSKIntegrals");
  if (func != NULL) {
    (*func)(sk, dist, sp1, sp2);
    return 1;
  }

  return 0;
}
