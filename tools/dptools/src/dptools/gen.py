#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2017  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

############################################################################
# Representation of the GEN format
############################################################################
import numpy as np
from dptools.common import openfile
from dptools.geometry import Geometry

__all__ = [ "Gen", ]

class Gen:
    """Representation of a GEN file.

    Attributes:
        geometry: Geometry object with atom positions and lattice vectors.
    """

    def __init__(self, geometry):
        """Initializes the instance.

        Args:
            geometry: geometry object containing the geometry information.
        """
        self.geometry = geometry


    @classmethod
    def fromfile(cls, fobj):
        """Creates a Gen instance from a file.

        Args:
            fobj: filename or file like object containing geometry in
                GEN-format.
        """
        fp = openfile(fobj, "r")
        lines = fp.readlines()
        fp.close()
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
        return cls(geometry)


    def tofile(self, fobj, relcoords=False):
        """Writes a GEN file.

        Args:
            fobj: File name or file object where geometry should be written.
            relcoords: If true, geometry will be outputted in relative
                coordinates (as multiples of the lattice vectors)
        """
        fp = openfile(fobj, "w")
        line = [ "{0:d}".format(self.geometry.natom), ]
        geo = self.geometry
        if geo.periodic:
            if relcoords:
                line.append("F")
                coords = geo.relcoords
            else:
                line.append("S")
                coords = geo.coords
        else:
            line.append("C")
            coords = geo.coords
        fp.write(" ".join(line) + "\n")
        fp.write(" ".join(geo.specienames) + "\n")
        for ii in range(geo.natom):
            fp.write("{0:6d} {1:3d} {2:18.10E} {3:18.10E} {4:18.10E}\n".format(
                ii + 1, geo.indexes[ii] + 1, *coords[ii]))
        if geo.periodic:
            fp.write("{0:18.10E} {1:18.10E} {2:18.10E}\n".format(*geo.origin))
            for vec in geo.latvecs:
                fp.write("{0:18.10E} {1:18.10E} {2:18.10E}\n".format(*vec))
        fp.close()
