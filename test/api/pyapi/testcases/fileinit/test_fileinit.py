#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

'''
Tests subsequent initializations/finalizations of the DFTB+ object by
calculating two incompatible geometries. Si2 and H2O geometries are
calculated alternately. After finishing a calculation, the finalization
routine of the interface will get called and a new geometry is set up
and calculated.
'''


import numpy as np
import dftbplus
from testhelpers import write_autotest_tag


LIB_PATH = '../../../../../prog/dftb+/libdftbplus'

# number of atoms (should match dftb_in.hsd files)
NATOM_SI = 2
NATOM_H2O = 3
MAX_ATOMS = 3

# number of iterations to calculate
NITER = 3


def init_collective_variables():
    '''Initialize collective variables, which will
       hold summed up results of multiple test runs.
    '''

    merminentot = 0.0
    gradientstot = np.zeros((MAX_ATOMS, 3))
    grosschgstot = np.zeros(MAX_ATOMS)

    return merminentot, gradientstot, grosschgstot

def update_collective_variables(merminen, gradients, grosschgs,
                                merminentot, gradientstot, grosschgstot):
    '''Add new results to collective variables.

    Args:
        merminen     (float):   mermin free energy
        gradients    (2darray): gradients
        grosschgs    (1darray): Gross charges
        merminentot  (float):   summation of mermin free energies
        gradientstot (2darray): summation of gradients
        grosschgstot (1darray): summation of Gross charges

    Returns:
        merminentot  (float):   updated mermin free energy
        gradientstot (2darray): updated gradients
        grosschgstot (1darray): updated Gross charges
    '''

    merminentot += merminen
    gradientstot = np.add(gradients, gradientstot)
    grosschgstot = np.add(grosschgs, grosschgstot)

    return merminentot, gradientstot, grosschgstot


def main():
    '''Main driver routine.'''

    # coordinates of Si2, in atomic units
    coords_si2 = np.array([
        [0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
        [2.2639291987021915, 2.4639291987021915, 2.5639291987021915]])

    # lattice vectors of Si2, in atomic units
    latvecs_si2 = np.array([
        [5.2278583974043830, 5.1278583974043830, 0.0000000000000000],
        [0.0000000000000000, 5.3278583974043830, 5.1278583974043830],
        [5.1278583974043830, 0.0000000000000000, 5.4278583974043830]])

    # coordinates of H2O, in atomic units
    coords_h2o = np.array([
        [0.00000000000E+00, -0.10000000000E+01,  0.00000000000E+00],
        [0.00000000000E+00,  0.00000000000E+00,  0.88306400000E+00],
        [0.00000000000E+00,  0.00000000000E+00, -0.78306400000E+00]])

    # initialize collective variables
    merminentot, gradientstot, grosschgstot = init_collective_variables()


    # dummy loop to test subsequent initializations/
    # finalizations of the DFTB+ object
    for iteration in range(0, NITER):

        tsi2 = not bool(iteration % 2)

        # use input for Si2 and H2O alternatingly, starting with Si2
        if tsi2:

            # set expected number of atoms
            natom0 = NATOM_SI

            cdftb = dftbplus.DftbPlus(libpath=LIB_PATH,
                                      hsdpath='dftb_in.Si2.hsd',
                                      logfile='log.Si2.log')

            # set geometry
            cdftb.set_geometry(coords_si2, latvecs=latvecs_si2)

        else:

            # set expected number of atoms
            natom0 = NATOM_H2O

            cdftb = dftbplus.DftbPlus(libpath=LIB_PATH,
                                      hsdpath='dftb_in.H2O.hsd',
                                      logfile='log.H2O.log')

            # set geometry
            cdftb.set_geometry(coords_h2o, latvecs=None)

        # get number of atoms
        natoms = cdftb.get_nr_atoms()

        # calculate energy, forces and Gross charges
        merminen = cdftb.get_energy()
        gradients = cdftb.get_gradients()
        grosschgs = cdftb.get_gross_charges()

        # finalize DFTB+ and clean up
        cdftb.close()

        if tsi2:

            dummygrads = np.zeros(3)
            dummychgs = 0
            merminentot, gradientstot, grosschgstot = \
            update_collective_variables(merminen,
                                        np.vstack((gradients, dummygrads)),
                                        np.hstack((grosschgs, dummychgs)),
                                        merminentot, gradientstot, grosschgstot)

        else:

            merminentot, gradientstot, grosschgstot = \
            update_collective_variables(merminen, gradients, grosschgs,
                                        merminentot, gradientstot, grosschgstot)

        if tsi2:

            # check whether calculator was initialized with correct nr. of atoms
            print('(Si2) Obtained nr. of atoms: {:d}'.format(natoms))
            print('(Si2) Expected nr. of atoms: {:d}\n'.format(natom0))

            # evaluate mermin free energy
            print('(Si2) Obtained Mermin-energy: ' + \
                  '{:15.10f}'.format(merminen))
            print('(Si2) Expected Mermin-energy: ' + \
                  '{:15.10f}\n'.format(-2.5897497363))

            # evaluate gradients
            print('(Si2) Obtained gradient of atom 1: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}'
                  .format(*gradients[0]))
            print('(Si2) Expected gradient of atom 1: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}'
                  .format(0.0306186399, 0.0026710677, -0.0007231241))
            print('(Si2) Obtained gradient of atom 2: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}'
                  .format(*gradients[1]))
            print('(Si2) Expected gradient of atom 2: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}\n'
                  .format(-0.0306186399, -0.0026710677, 0.0007231241))

            # evaluate Gross charges
            print('(Si2) Obtained Gross charges: {:15.10f} {:15.10f}'
                  .format(*grosschgs))
            print('(Si2) Expected Gross charges: {:15.10f} {:15.10f}\n\n'
                  .format(0.0, 0.0))

        else:

            # check whether calculator was initialized with correct nr. of atoms
            print('(H2O) Obtained nr. of atoms: {:d}'.format(natoms))
            print('(H2O) Expected nr. of atoms: {:d}\n'.format(natom0))

            # evaluate mermin free energy
            print('(H2O) Obtained Mermin-energy: ' + \
                  '{:15.10f}'.format(merminen))
            print('(H2O) Expected Mermin-energy: ' + \
                  '{:15.10f}\n'.format(-3.7534584715))

            # evaluate gradients
            print('(H2O) Obtained gradient of atom 1: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}'
                  .format(*gradients[0]))
            print('(H2O) Expected gradient of atom 1: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}'
                  .format(0.0178274222, 1.1052225107, -0.0835776976))
            print('(H2O) Obtained gradient of atom 2: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}'
                  .format(*gradients[1]))
            print('(H2O) Expected gradient of atom 2: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}'
                  .format(-0.0073820765, -0.4646991011, -0.5429357096))
            print('(H2O) Obtained gradient of atom 3: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}'
                  .format(*gradients[2]))
            print('(H2O) Expected gradient of atom 3: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}\n'
                  .format(-0.0058552790, -0.6382027486, 0.6279471931))

            # evaluate Gross charges
            print('(H2O) Obtained Gross charges: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}'
                  .format(*grosschgs))
            print('(H2O) Expected Gross charges: ' + \
                  '{:15.10f} {:15.10f} {:15.10f}\n\n'
                  .format(-0.6519945363, 0.3314953102, 0.3204992261))


    # --------------------------WRITE AUTOTEST.TAG------------------------------

    # write autotest.tag file, containing the collective variables
    write_autotest_tag('autotest.tag', freeEgy=merminentot,
                       forceTot=-gradientstot, qOutAtGross=grosschgstot)


if __name__ == "__main__":
    main()
