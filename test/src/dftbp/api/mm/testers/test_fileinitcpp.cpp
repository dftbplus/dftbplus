/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2025  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

#include <iostream>
#include "dftbplus.h"
#include "testhelpers.h"
using namespace std;

// Reads in a dftb_in.hsd file before calculating properties.

int main() {
  DftbPlus calculator;
  DftbPlusInput input;

  double mermin_energy;
  int major, minor, patch;
  string filename;

  dftbp_api(&major, &minor, &patch);
  cout << "DFTB+ API version " << major << "." << minor << "." << patch << "\n";

  bool instsafe = dftbp_is_instance_safe();
  printf(instsafe ? "API is instance safe\n" : "API is NOT instance safe\n");

  // Initialize a DFTB+ calculator
  dftbp_init(&calculator, NULL);

  cout << "Input file name?\n";
  cin >> filename;

  // Initialize DFTB+ input tree from input in external file
  dftbp_get_input_from_file(&calculator, filename.c_str(), &input);

  dftbp_process_input(&calculator, &input);
  dftbp_input_final(&input);

  int natom = dftbp_get_nr_atoms(&calculator);
  cout << "Obtained nr. of atoms: " << natom << "\n";

  dftbp_get_energy(&calculator, &mermin_energy);

  cout.precision(8);

  cout << "Mermin free energy: " << mermin_energy << " H\n";

  dftbp_write_autotest_tag(natom, 0, 0, mermin_energy, NULL, NULL, NULL, NULL, NULL, NULL);

  return 0;
}
