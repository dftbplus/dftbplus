#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

'''
Tests the scenario that a geometry gets slightly displaced after initialization.
'''


import numpy as np
import dftbplus
from testhelpers import write_autotest_tag


LIB_PATH = '../../../../../prog/dftb+/libdftbplus'

# number of atoms (should match dftb_in.hsd)
NATOM0 = 2

# DFTB+ conversion factors
# (according to prog/dftb+/lib_common/constants.F90)
BOHR__AA = 0.529177249
AA__BOHR = 1 / BOHR__AA


def main():
    '''Main driver routine.'''

    # initial coordinates of Si2, in atomic units
    initialcoords = np.array([
        [0.0000000000000000, 0.0000000000000000, 0.0000000000000000],
        [2.5639291987021915, 2.5639291987021915, 2.5639291987021915]])

    # small displacement in the coordinates
    coords = np.copy(initialcoords)
    coords[0, 0] = initialcoords[0, 0] + 0.2 * AA__BOHR

    # initial lattice vectors of Si2, in atomic units
    initiallatvecs = np.array([
        [5.1278583974043830, 5.1278583974043830, 0.0000000000000000],
        [0.0000000000000000, 5.1278583974043830, 5.1278583974043830],
        [5.1278583974043830, 0.0000000000000000, 5.1278583974043830]])

    # small displacement in the lattice vectors
    latvecs = np.copy(initiallatvecs)
    latvecs[0, 0] = initiallatvecs[0, 0] + 0.1 * AA__BOHR


    # ------------------CALCULATION OF INITIAL GEOMETRY-------------------------

    cdftb = dftbplus.DftbPlus(libpath=LIB_PATH,
                              hsdpath='dftb_in.hsd',
                              logfile='log.log')

    # set initial geometry
    cdftb.set_geometry(initialcoords, latvecs=initiallatvecs)

    # get number of atoms
    natoms = cdftb.get_nr_atoms()

    # calculate energy, forces and Gross charges
    merminen = cdftb.get_energy()
    gradients = cdftb.get_gradients()
    grosschg = cdftb.get_gross_charges()

    # check whether calculator was initialized with correct nr. of atoms
    print('(Si2) Obtained nr. of atoms: {:d}'.format(natoms))
    print('(Si2) Expected nr. of atoms: {:d}\n'.format(NATOM0))

    # evaluate mermin free energy
    print('(Si2) Obtained Mermin-energy: {:15.10f}'.format(merminen))
    print('(Si2) Expected Mermin-energy: {:15.10f}\n'.format(-2.5933460731))

    # evaluate gradients
    print('(Si2) Obtained gradient of atom 1: {:15.10f} {:15.10f} {:15.10f}'
          .format(*gradients[0]))
    print('(Si2) Expected gradient of atom 1: {:15.10f} {:15.10f} {:15.10f}'
          .format(-0.0103215090, -0.0103215090, -0.0103215090))
    print('(Si2) Obtained gradient of atom 2: {:15.10f} {:15.10f} {:15.10f}'
          .format(*gradients[1]))
    print('(Si2) Expected gradient of atom 2: {:15.10f} {:15.10f} {:15.10f}\n'
          .format(0.0103215090, 0.0103215090, 0.0103215090))

    # evaluate Gross charges
    print('(Si2) Obtained Gross charges: {:15.10f} {:15.10f}'
          .format(*grosschg))
    print('(Si2) Expected Gross charges: {:15.10f} {:15.10f}\n\n'
          .format(0.0, 0.0))


    # ------------------CALCULATION OF DISPLACED GEOMETRY-----------------------

    # set displaced geometry
    cdftb.set_geometry(coords, latvecs=latvecs)

    # get number of atoms
    natoms = cdftb.get_nr_atoms()

    # calculate energy, forces and Gross charges
    merminen = cdftb.get_energy()
    gradients = cdftb.get_gradients()
    grosschg = cdftb.get_gross_charges()

    # check whether calculator was initialized with correct nr. of atoms
    print('(Si2) Obtained nr. of atoms: {:d}'.format(natoms))
    print('(Si2) Expected nr. of atoms: {:d}\n'.format(NATOM0))

    # evaluate mermin free energy
    print('(Si2) Obtained Mermin-energy: {:15.10f}'.format(merminen))
    print('(Si2) Expected Mermin-energy: {:15.10f}\n'.format(-2.5854559427))

    # evaluate gradients
    print('(Si2) Obtained gradient of atom 1: {:15.10f} {:15.10f} {:15.10f}'
          .format(*gradients[0]))
    print('(Si2) Expected gradient of atom 1: {:15.10f} {:15.10f} {:15.10f}'
          .format(0.0442660049, -0.0147463633, -0.0193148538))
    print('(Si2) Obtained gradient of atom 2: {:15.10f} {:15.10f} {:15.10f}'
          .format(*gradients[1]))
    print('(Si2) Expected gradient of atom 2: {:15.10f} {:15.10f} {:15.10f}\n'
          .format(-0.0442660049, 0.0147463633, 0.0193148538))

    # evaluate Gross charges
    print('(Si2) Obtained Gross charges: {:15.10f} {:15.10f}'
          .format(*grosschg))
    print('(Si2) Expected Gross charges: {:15.10f} {:15.10f}'
          .format(0.0, 0.0))

    # finalize DFTB+ and clean up
    cdftb.close()


    # --------------------------WRITE AUTOTEST.TAG------------------------------

    # write autotest.tag file of water molecule calculation
    write_autotest_tag('autotest.tag', freeEgy=merminen,
                       forceTot=-gradients, qOutAtGross=grosschg)


if __name__ == "__main__":
    main()
