#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''Information in CIF-files (only basic informations)'''

import numpy as np
from dptools.common import openfile
from dptools.geometry import Geometry
from dptools.geometry import get_latvecs_fromcif

__all__ = ['Cif']


_ABSTOLERANCE = 1E-10
_RELTOLERANCE = 1E-10


class Cif:
    '''Representation of a CIF file.

    Attributes:
        geometry: Geometry object with atom positions and lattice vectors.
        celllengths: Length of the 3 cell vectors.
        cellangles: Angles between the cell vectors in radian.
    '''

    def __init__(self, geometry):
        '''Initializes a CIF instance.

        Args:
            geometry: geometry object with atom positions and lattice vectors.
        '''
        self.geometry = geometry
        self.celllengths = np.array([np.sqrt(np.sum(vv**2))
                                     for vv in geometry.latvecs], dtype=float)
        # cellangles in radians (alpha, beta and gamma as in crystallography)
        self.cellangles = np.empty(3, dtype=float)
        for ii in range(3):
            i1 = (ii + 1) % 3
            i2 = (ii + 2) % 3
            v1 = geometry.latvecs[i1]
            v2 = geometry.latvecs[i2]
            dot = np.dot(v1, v2) / (self.celllengths[i1] * self.celllengths[i2])
            self.cellangles[ii] = np.arccos(dot)

    @classmethod
    def fromfile(cls, fobj):
        '''Reads crystallographic information from a CIF file.

        Args:
            fobj: filename or file like object containing geometry in
                CIF-format.

        '''
        fp = openfile(fobj, "r")
        lines = fp.readlines()
        fp.close()
        celllengths = np.empty(3, dtype=float)
        cellangles = np.empty(3, dtype=float)
        for jj in range(3):
            celllengths[jj] = lines[jj + 1].split()[1]
            cellangles[jj] = lines[jj + 4].split()[1]
        natom = len(lines) - 13
        specienames = []
        speciedict = {}
        indexes = np.empty((natom, ), dtype=int)
        coords = np.empty((natom, 3), dtype=float)
        for ii, line in enumerate(lines[13:13+natom]):
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
        latvecs = get_latvecs_fromcif(celllengths, cellangles)
        geometry = Geometry(specienames, indexes, coords, latvecs=latvecs, relcoords=True)
        return cls(geometry)


    def tofile(self, fobj):
        '''Writes a CIF file.

        Args:
            fobj: File name or file object where geometry should be written.
        '''
        geo = self.geometry
        fp = openfile(fobj, "w")
        fp.write("data_global\n")
        for name, value in zip(["a", "b", "c"], self.celllengths):
            fp.write("_cell_length_{0:s} {1:.10f}\n".format(name, value))
        # cell angles are needed in degrees
        for name, value in zip(["alpha", "beta", "gamma"],
                               self.cellangles * 180.0 / np.pi):
            fp.write("_cell_angle_{0:s} {1:.10f}\n".format(name, value))
        fp.write("_symmetry_space_group_name_H-M 'P 1'\n")
        fp.write("loop_\n_atom_site_label\n_atom_site_fract_x\n"
                 "_atom_site_fract_y\n_atom_site_fract_z\n")
        for ii in range(geo.natom):
            fp.write("{0:3s} {1:.10f} {2:.10f} {3:.10f}\n".format(
                geo.specienames[geo.indexes[ii]], *geo.relcoords[ii]))
        fp.close()

    def equals(self, other, abstolerance=_ABSTOLERANCE, reltolerance=_RELTOLERANCE):
        '''Checks whether object equals to an other one.

        Args:
            other (Cif): Other Cif object.
            tolerance (float): Maximal allowed deviation in floating point
                numbers (e.g. coordinates).
        '''
        celllengths_close = np.allclose(self.celllengths, other.celllengths,
                                        rtol=reltolerance, atol=abstolerance)
        cellangles_close = np.allclose(self.cellangles, other.cellangles,
                                       rtol=reltolerance, atol=abstolerance)
        if not celllengths_close:
            return False
        if not cellangles_close:
            return False
        if not self.geometry.equals(other.geometry, abstolerance):
            return False
        return True
