/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2023  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

/**
 * Explicitly sets the neighbour list via the API and checks the result
 */

#include <math.h>
#include <stdio.h>
#include "dftbplus.h"
#include "testhelpers.h"

int main()
{
  DftbPlus calculator;
  DftbPlusInput input;

  int major, minor, patch;
  dftbp_api(&major, &minor, &patch);
  printf("API version %d.%d.%d\n", major, minor, patch);

  dftbp_init(&calculator, NULL);

  dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);
  dftbp_process_input(&calculator, &input);

  // setup all data for the neighbour list
  double cutoff = 6.0;
  double x = 2.5639291987021915;
  int nNeighbour[2] = {4.0, 0.0};
  int iNeighbour[8] = {1, 2, 3, 4, 0, 0, 0, 0};
  double dist = sqrt(19.721198807872987);
  double neighDist[8] = {dist, dist, dist, dist, 0.0, 0.0, 0.0, 0.0};
  double coord[12] = {-x, -x, x, x, -x, -x, -x, x, -x, x, x, x};
  int img2CentCell[4] = {2, 2, 2, 2};

  dftbp_set_neighbour_list(&calculator, 4, 4, nNeighbour, iNeighbour, neighDist, cutoff,
                           coord, img2CentCell);

  // evaluate energy and forces
  double mermin_energy;
  dftbp_get_energy(&calculator, &mermin_energy);

  double gradients[6];
  dftbp_get_gradients(&calculator, gradients);

  dftbp_final(&calculator);

  dftbp_write_autotest_tag(2, 0, 0, mermin_energy, gradients, NULL, NULL, NULL, NULL, NULL);

  return 0;
}
