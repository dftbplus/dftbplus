/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2024  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/
#ifndef __DFTBPLUS_PLUGIN_H__
#define __DFTBPLUS_PLUGIN_H__

#ifdef __cplusplus
extern "C" {
#endif

void init();
void final();

void getSKIntegrals(double*, double, int, int, int, int);

#ifdef __cplusplus
}
#endif

#endif
