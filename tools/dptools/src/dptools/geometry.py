#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Representation of a geometry'''

import numpy as np
import numpy.linalg as la

__all__ = ["Geometry"]


class Geometry:
    """Atomic geometry representation.

    Attributes:
        specienames: Name of atomtypes which can be found in the geometry.
        nspecie: Number of species.
        indexes: For each atom the index of the corresponding specie
            in specienames.
        natom: Number of atoms.
        coords: xyz coordinates of the atoms.
        origin: Origin.
        periodic: True if the structure is periodic.
        latvecs: Lattice vectors (None for non-periodic structures).
        relcoords: Relative lattice coordinates (None if non-periodic).
    """

    def __init__(self, specienames, indexes, coords, latvecs=None, origin=None,
                 relcoords=False):
        """Initializes a geometry object.

        Args:
            specienames: Names of the species occuring in the geometry.
            indexes: Species index for every atom. Shape (natom,)
            coords: Coordinates of the atoms.
            latvecs: Lattice vectors (default: None, non-periodic modell).
            origin: Origin of the primitive cell (default (0.0, 0.0, 0.0))
            relcoords: If set to yes, coordinates are assumed to be relative
                coordinates specified as multiples of the lattice vectors.
        """
        self.specienames = list(specienames)
        self.nspecie = len(self.specienames)
        self.indexes = np.array(indexes)
        self.natom = len(self.indexes)
        self.coords = np.array(coords)
        if origin is None:
            self.origin = np.array((0, 0, 0), dtype=float)
        else:
            self.origin = np.array(origin, dtype=float)
        self.periodic = latvecs is not None
        if self.periodic:
            self.latvecs = np.array(latvecs, dtype=float)
            self._invlatvecs = la.inv(self.latvecs)
            if relcoords:
                self.relcoords = coords
                self.coords = np.dot(self.relcoords, self.latvecs) + self.origin
            else:
                self.coords = coords
                self.relcoords = np.dot(self.coords - self.origin,
                                        self._invlatvecs)
        else:
            self.latvecs = None
            self._invlatvecs = None
            self.relcoords = None


    def setlattice(self, latvecs, origin=None):
        """Makes geometry periodic or changes supercell vectors.

        Args:
            latvecs: Periodicity defined by lattice vectors.
            origin: Origin (default: 0, 0, 0)
        """
        self.latvecs = np.array(latvecs, dtype=float)
        self._invlatvecs = la.inv(self.latvecs)
        if origin is None:
            self.origin = np.array((0, 0, 0), dtype=float)
        else:
            self.origin = np.array(origin, dtype=float)
        self.relcoords = np.dot(self.coords - self.origin, self._invlatvecs)
        self.periodic = True


    def equals(self, other, tolerance):
        '''Checks whether object equals to an other one.

        Args:
            other (Geometry): Other geometry.
            tolerance (float): Maximal allowed deviation in floating point
                numbers (e.g. coordinates).
        '''
        if self.specienames != other.specienames:
            return False
        if np.any(self.indexes != other.indexes):
            return False
        if np.any(abs(self.coords - other.coords) > tolerance):
            return False
        if np.any(abs(self.origin - other.origin) > tolerance):
            return False
        if self.periodic != other.periodic:
            return False
        if self.periodic:
            if np.any(abs(self.latvecs - other.latvecs) > tolerance):
                return False
            if np.any(abs(self.relcoords - other.relcoords) > tolerance):
                return False
        return True


def get_latvecs_fromcif(celllengths, cellangles):
    '''Calculate cartesian lattice vectors from crystallographic CIF information

    Args:
        celllengths: Cell lengths a, b, c from CIF
        cellangles: Angles alpha, beta, gamma that span the cell

    Returns:
        latvecs: calculated cartesian lattice vectors
    '''
    celllengths = celllengths
    cellangles = cellangles * 2*np.pi/360
    latvecs = np.empty((3, 3), dtype=float)

    latvecs[0, :] = np.array([celllengths[0], 0, 0], dtype=float)
    latvecs[1, :] = np.array([np.cos(cellangles[2]), np.sin(cellangles[2]), 0], dtype=float) \
    * celllengths[1]
    trig_term = (np.cos(cellangles[0]) - np.cos(cellangles[2]) * np.cos(cellangles[1])) \
    /np.sin(cellangles[2])
    latvecs[2, :] = np.array([np.cos(cellangles[1]), trig_term, \
    np.sqrt(1 - np.cos(cellangles[1])**2 - trig_term**2)], dtype=float) * celllengths[2]

    return latvecs
