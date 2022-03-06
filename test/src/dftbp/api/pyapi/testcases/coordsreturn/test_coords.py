#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

'''
Tests return of coordinates from the DFTB+ object and their use in
evaluating some properties by calculating two incompatible
geometries. Si2 and H2O geometries are calculated alternately.
'''

import numpy as np
import dftbplus
from testhelpers import write_autotest_tag

LIB_PATH = '../../../../../../../src/dftbp/libdftbplus'

# Finite difference step
delta = 0.000005

def main():
    '''Main driver routine.'''

    latvecs_si2 = np.array([
        [5.2278583974043830, 5.1278583974043830, 0.0000000000000000],
        [0.0000000000000000, 5.3278583974043830, 5.1278583974043830],
        [5.1278583974043830, 0.0000000000000000, 5.4278583974043830]])

    inputcase = 'H2O'
    dftbcalc = dftbplus.DftbPlus(libpath=LIB_PATH,
                                 hsdpath='dftb_in.'+inputcase+'.hsd',
                                 logfile=inputcase+'.log')

    # get number of atoms and coordinates
    natoms = dftbcalc.get_nr_atoms()
    coords = dftbcalc.get_coords()

    mu = np.zeros(3)
    zstar = np.zeros((3, natoms, 3))
    qdot = np.zeros((3, natoms))

    for icart in range(3):
        energy = np.zeros(2)
        gradients = np.zeros((2, natoms, 3))
        grosschgs = np.zeros((2, natoms))
        for jj in range(2):
            efield = np.zeros(3)
            deltaE = float(2 * jj - 1) * delta
            efield[icart] = deltaE

            extpot = np.dot(efield, coords.T)
            extpotgrad = np.zeros((natoms, 3))
            extpotgrad[:,icart] = deltaE

            # set external potential
            dftbcalc.set_external_potential(extpot, extpotgrad=extpotgrad)

            # set geometry, as work-around for API bug
            if inputcase == 'Si2':
                dftbcalc.set_geometry(coords, latvecs_si2)
            else:
                dftbcalc.set_geometry(coords, None)

            # calculate energy, forces and gross (Mulliken) charges
            energy[jj] = dftbcalc.get_energy()
            gradients[jj] = dftbcalc.get_gradients()
            grosschgs[jj] = dftbcalc.get_gross_charges()

        mu[icart] = (energy[0] - energy[1]) / (2.0 * delta)
        zstar[icart] = (gradients[0] - gradients[1]) / (2.0 * delta)
        qdot[icart] = (grosschgs[0] - grosschgs[1]) / (2.0 * delta)

    # finalize DFTB+ and clean up
    dftbcalc.close()

    print('dipole',mu)
    print('Born')
    for iat in range(natoms):
        print("atom %i" % (iat+1))
        for icart in range(3):
            print(zstar[icart][iat][:])
    print('Atom polarisability',qdot.T)

if __name__ == "__main__":
    main()
