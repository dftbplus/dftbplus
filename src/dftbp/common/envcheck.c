/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2023  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

#include <sys/resource.h>

/**
 * Queries current stacksize.
 *
 * \param[out] cstack Current stacksize of task.
 *
 * \param[out] ierr Error status.
 */
void get_stacksize_c(int *cstack, int *ierr) {

  struct rlimit rlim;
  *ierr = getrlimit(RLIMIT_STACK, &rlim);

  *cstack = rlim.rlim_cur;

}
