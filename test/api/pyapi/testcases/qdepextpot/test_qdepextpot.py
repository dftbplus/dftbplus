#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

'''
Tests the capability of the interface to run a calculation with a population
dependent external potential. Since only the user knows how to calculate it,
callback functions can be defined, which will be executed at runtime.
'''


import numpy as np
import dftbplus
from testhelpers import write_autotest_tag


LIB_PATH = '../../../../../prog/dftb+/libdftbplus'

# number of atoms (should match dftb_in.hsd file)
NATOM = 3


class PotentialCalculator:
    '''

       Auxiliary class for calculating the population dependent external
       potential and its gradients. An instance of this class gets handed over
       to DFTB+ via the ctypes interface, to handle the necessary callbacks.

    '''


    def __init__(self, qmcoords, mmcoords, mmcharges):
        '''Initializes a PotentialCalculator object.

        Args:

            qmcoords  (2darray): coordinates of QM-atoms (shape: [qmatoms, 3])
            mmcoords  (2darray): coordinates of MM-atoms (shape: [mmatoms, 3])
            mmcharges (1darray): charges of MM-atoms     (shape: [mmatoms, 1])

        '''

        self._qmcoords = qmcoords
        self._mmcoords = mmcoords

        self._qmatoms = np.shape(self._qmcoords)[0]
        self._mmatoms = np.shape(self._mmcoords)[0]

        self._mmcharges = mmcharges


    def calc_extpot(self, dqatom):
        '''Calculates the current external potential
           using the properties of the MM- and QM-atoms.

        Args:

            dqatom (1darray): population difference with respect to reference
                              population (usually the neutral atom)
                              Note: population means electrons, so a positive
                              number indicates electron excess

        Returns:

            extpot (1darray): updated external potential at
                              the position of each QM-atom

        '''

        # Note: Some types of potential require knowledge of the
        # current atomic populations, which is provided by dqatom.

        extpot = np.zeros(self._qmatoms)

        for iatqm in range(0, self._qmatoms):
            extpot[iatqm] = 0.0
            qmpos = self._qmcoords[iatqm, :]
            for iatmm in range(0, self._mmatoms):
                mmpos = self._mmcoords[iatmm, :]
                mmcharge = self._mmcharges[iatmm]
                dist = distance(qmpos, mmpos)
                extpot[iatqm] += -mmcharge / dist

        return extpot


    def calc_extpotgrad(self, dqatom):
        '''Calculates the current gradients of the external
           potential using the properties of the MM- and QM-atoms.

        Args:

            dqatom     (1darray): population difference with respect
                                  to reference population
                                  (usually the neutral atom)
                                  Note: population means electrons, so a
                                  positive number indicates electron excess

        Returns:

            extpotgrad (2darray): updated potential gradient at
                                  the position of each QM-atom

        '''

        # Note: Some types of potential require knowledge of the
        # current atomic populations, which is provided by dqatom.

        extpotgrad = np.zeros((self._qmatoms, 3))

        for iatqm in range(0, self._qmatoms):
            qmpos = self._qmcoords[iatqm, :]
            extpotgrad[iatqm, :] = 0.0
            for iatmm in range(0, self._mmatoms):
                mmpos = self._mmcoords[iatmm, :]
                mmcharge = self._mmcharges[iatmm]
                dist = distance(qmpos, mmpos)
                dist3 = dist**3
                extpotgrad[iatqm, 0] += -mmcharge * (mmpos[0]-qmpos[0]) / dist3
                extpotgrad[iatqm, 1] += -mmcharge * (mmpos[1]-qmpos[1]) / dist3
                extpotgrad[iatqm, 2] += -mmcharge * (mmpos[2]-qmpos[2]) / dist3

        return extpotgrad


def get_extpot(potcalc, dqatom, extpotatom):
    '''Queries the external potential.

    Args:

        potcalc   (pyobject): instance of a class that provides
                              methods for calculating the external
                              potential and its gradients
        dqatom     (1darray): population difference with respect to reference
                              population (usually the neutral atom)
                              Note: population means electrons, so a positive
                              number indicates electron excess
        extpotatom (1darray): potential at the position of each QM-atom
                              Note: it should be the potential as felt by
                              an electron (negative potential value means
                              attraction for an electron)

    '''

    extpotatom[:] = potcalc.calc_extpot(dqatom)


def get_extpotgrad(potcalc, dqatom, extpotatomgrad):
    '''Queries the external potentials gradients.

    Args:

        potcalc       (pyobject): instance of a class that provides
                                  methods for calculating the external
                                  potential and its gradients
        dqatom         (1darray): population difference with respect to
                                  reference population
                                  (usually the neutral atom)
                                  Note: population means electrons, so a
                                  positive number indicates electron excess
        extpotatomgrad (2darray): potential gradient at the position of each
                                  QM-atom
                                  Note: it should be the gradient of the
                                  potential as felt by an electron (negative
                                  potential value means attraction for an
                                  electron)
    '''

    extpotatomgrad[:, :] = potcalc.calc_extpotgrad(dqatom)


def distance(aa, bb):
    '''Auxiliary function for calculating
       the distance between two vectors.

    Args:

       aa (1darray): vector a
       bb (1darray): vector b

    Returns:

       dist (float): distance between vectors

    '''

    aa = np.asarray(aa)
    bb = np.asarray(bb)

    dist = np.linalg.norm(aa - bb)

    return dist


def main():
    '''Main driver routine.'''

    # coordinates of H2O, in atomic units
    qmcoords = np.array([
        [0.0000000000000000, -1.8897259885789233,  0.0000000000000000],
        [0.0000000000000000,  0.0000000000000000,  1.4797763915205659],
        [0.0000000000000000,  0.0000000000000000, -1.4797763915205659]])

    # coordinates of MM-charges, in atomic units
    mmcoords = np.array([
        [-0.944863438887178, -9.44863438887178, 1.70075418999692],
        [ 4.34637181888102,  -5.85815332110050, 2.64561762888410]])

    # MM-charges, in atomic units
    mmcharges = np.array([2.5, -1.9])

    potcalc = PotentialCalculator(qmcoords, mmcoords, mmcharges)


    # ----------------------------MAIN-CALCULATION------------------------------

    cdftb = dftbplus.DftbPlus(libpath=LIB_PATH,
                              hsdpath='dftb_in.hsd',
                              logfile='log.log')

    # set geometry
    cdftb.set_geometry(qmcoords, latvecs=None)

    cdftb.register_ext_pot_generator(potcalc, get_extpot, get_extpotgrad)

    # get number of atoms
    natoms = cdftb.get_nr_atoms()

    # calculate energy, forces and Gross charges
    merminen = cdftb.get_energy()
    gradients = cdftb.get_gradients()
    grosschgs = cdftb.get_gross_charges()

    # check whether calculator was initialized with correct nr. of atoms
    print('(H2O) Obtained nr. of atoms: {:d}'.format(natoms))
    print('(H2O) Expected nr. of atoms: {:d}\n'.format(NATOM))

    # evaluate mermin free energy
    print('(H2O) Obtained Mermin-energy: ' +
          '{:15.10f}'.format(merminen))
    print('(H2O) Expected Mermin-energy: ' +
          '{:15.10f}\n'.format(-3.9854803392))

    # evaluate gradients
    print('(H2O) Obtained gradient of atom 1: ' +
          '{:15.10f} {:15.10f} {:15.10f}'
          .format(*gradients[0]))
    print('(H2O) Expected gradient of atom 1: ' +
          '{:15.10f} {:15.10f} {:15.10f}'
          .format(0.0176513638, -0.1831376018, 0.0031982515))
    print('(H2O) Obtained gradient of atom 2: ' +
          '{:15.10f} {:15.10f} {:15.10f}'
          .format(*gradients[1]))
    print('(H2O) Expected gradient of atom 2: ' +
          '{:15.10f} {:15.10f} {:15.10f}'
          .format(-0.0061402266, 0.0955090293, 0.0394035230))
    print('(H2O) Obtained gradient of atom 3: ' +
          '{:15.10f} {:15.10f} {:15.10f}'
          .format(*gradients[2]))
    print('(H2O) Expected gradient of atom 3: ' +
          '{:15.10f} {:15.10f} {:15.10f}\n'
          .format(-0.0037720260, 0.0923535862, -0.0402979580))

    # evaluate Gross charges
    print('(H2O) Obtained Gross charges: ' +
          '{:15.10f} {:15.10f} {:15.10f}'
          .format(*grosschgs))
    print('(H2O) Expected Gross charges: ' +
          '{:15.10f} {:15.10f} {:15.10f}\n\n'
          .format(-0.4943983279, 0.2641722128, 0.2302261151))

    # finalize DFTB+ and clean up
    cdftb.close()


    # --------------------------WRITE AUTOTEST.TAG------------------------------

    # write autotest.tag file of water molecule calculation
    write_autotest_tag('autotest.tag', freeEgy=merminen,
                       forceTot=-gradients, qOutAtGross=grosschgs)


if __name__ == "__main__":
    main()
