#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''Representation of the XYZ-format'''

import numpy as np
from dptools.common import openfile
from dptools.geometry import Geometry

__all__ = ["Xyz"]


_TOLERANCE = 1E-10


class Xyz:
    """Representation of an XYZ-file.

    Attributes:
        geometry: Geometry object with atom positions.
        comment: Content of the comment line in the XYZ-file.
    """

    def __init__(self, geometry, comment=""):
        """Creates XYZ-instance.

        Args:
            geometry: geometry object with atom positions.
            comment: additional comment (should be max. 1 line)
        """
        self.geometry = geometry
        self.comment = comment


    @classmethod
    def fromfile(cls, fobj):
        """Creates an XYZ instance from a file.

        Args:
            fobj: filename or file like object containing geometry in
                XYZ-format.

        """
        fp = openfile(fobj, "r")
        lines = fp.readlines()
        fp.close()
        words = lines[0].split()
        natom = int(words[0])
        comment = lines[1].strip()
        specienames = []
        speciedict = {}
        indexes = np.empty((natom, ), dtype=int)
        coords = np.empty((natom, 3), dtype=float)
        for ii, line in enumerate(lines[2:2+natom]):
            words = line.split()
            species = words[0]
            index = speciedict.get(species, -1)
            if index == -1:
                specienames.append(species)
                speciedict[species] = len(specienames) - 1
                indexes[ii] = len(specienames) - 1
            else:
                indexes[ii] = index
            coords[ii] = np.array(words[1:4], dtype=float)
        geometry = Geometry(specienames, indexes, coords)
        return cls(geometry, comment)


    def tofile(self, fobj):
        """Writes an XYZ file.

        Args:
            fobj: File name or file object where geometry should be written.
        """
        fp = openfile(fobj, "w")
        geo = self.geometry
        fp.write("{0:d}\n".format(geo.natom))
        fp.write(self.comment + "\n")
        for ii in range(geo.natom):
            fp.write("{0:3s} {1:18.10E} {2:18.10E} {3:18.10E}\n".format(
                geo.specienames[geo.indexes[ii]], *geo.coords[ii]))
        fp.close()

    def equals(self, other, tolerance=_TOLERANCE, check_comment=False):
        '''Checks whether object equals to an other one.

        Args:
            other (Xyz): Other Xyz object.
            tolerance (float): Maximal allowed deviation in floating point
                numbers (e.g. coordinates).
        '''
        if not self.geometry.equals(other.geometry, tolerance):
            return False
        if check_comment:
            if self.comment != other.comment:
                return False
        return True
