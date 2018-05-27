#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2018  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Representation of the GEN format.'''

import numpy as np
from dptools.common import OpenFile
from dptools.geometry import Geometry

__all__ = ["Gen"]


_TOLERANCE = 1E-10


class Gen:
    """Representation of a GEN file.

    Attributes:
        geometry: Geometry object with atom positions and lattice vectors.
    """

    def __init__(self, geometry, fractional=False):
        """Initializes the instance.

        Args:
            geometry: geometry object containing the geometry information.
            fractional: Whether fractional coordinates are preferred.
        """
        self.geometry = geometry
        self.fractional = fractional


    @classmethod
    def fromfile(cls, fobj):
        """Creates a Gen instance from a file.

        Args:
            fobj: filename or file like object containing geometry in
                GEN-format.
        """
        with OpenFile(fobj, 'r') as fp:
            lines = fp.readlines()
        words = lines[0].split()
        natom = int(words[0])
        flag = words[1].lower()
        if flag == "s":
            periodic = True
            relative = False
        elif flag == "f":
            periodic = True
            relative = True
        else:
            periodic = False
            relative = False
        specienames = lines[1].split()
        indexes = np.empty((natom, ), dtype=int)
        coords = np.empty((natom, 3), dtype=float)
        for ii, line in enumerate(lines[2:2+natom]):
            words = line.split()
            indexes[ii] = int(words[1]) - 1
            coords[ii] = np.array(words[2:5], dtype=float)
        if periodic:
            origin = np.array(lines[natom+2].split(), dtype=float)
            latvecs = np.empty((3, 3), dtype=float)
            for jj in range(3):
                latvecs[jj] = np.array(lines[natom+3+jj].split(), dtype=float)
        else:
            origin = None
            latvecs = None
        geometry = Geometry(specienames, indexes, coords, latvecs, origin,
                            relative)
        return cls(geometry, relative)


    def tofile(self, fobj):
        """Writes a GEN file.

        Args:
            fobj: File name or file object where geometry should be written.
        """
        lines = []
        line = ["{0:d}".format(self.geometry.natom)]
        geo = self.geometry
        if geo.periodic:
            if self.fractional:
                line.append("F")
                coords = geo.relcoords
            else:
                line.append("S")
                coords = geo.coords
        else:
            line.append("C")
            coords = geo.coords

        coords = _round_to_zero(coords, _TOLERANCE)
        lines.append(" ".join(line) + "\n")
        lines.append(" ".join(geo.specienames) + "\n")
        for ii in range(geo.natom):
            lines.append("{0:6d} {1:3d} {2:18.10E} {3:18.10E} {4:18.10E}\n"\
                         .format(ii + 1, geo.indexes[ii] + 1, *coords[ii]))
        if geo.periodic:
            origin = _round_to_zero(geo.origin, _TOLERANCE)
            lines.append("{0:18.10E} {1:18.10E} {2:18.10E}\n".format(*origin))
            latvecs = _round_to_zero(geo.latvecs, _TOLERANCE)
            for vec in latvecs:
                lines.append("{0:18.10E} {1:18.10E} {2:18.10E}\n".format(*vec))
        with OpenFile(fobj, 'w') as fp:
            fp.writelines(lines)


    def equals(self, other, tolerance=_TOLERANCE):
        '''Checks whether object equals to an other one.

        Args:
            other (Gen): Other Gen object.
            tolerance (float): Maximal allowed deviation in floating point
                numbers (e.g. coordinates).
        '''
        if not self.fractional == other.fractional:
            return False
        if not self.geometry.equals(other.geometry, tolerance):
            return False
        return True


def _round_to_zero(array, tolerance):
    '''Rounds elements of an array to zero below given tolerance.'''
    if tolerance is None:
        return array
    else:
        return np.where(abs(array) < tolerance, 0.0, array)
