#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''Representation of the NXY-format'''

import numpy as np

__all__ = ['Nxy']


_ABSTOLERANCE = 1E-10
_RELTOLERANCE = 1E-10


class Nxy:
    '''Representation of an NXY-file.

    Attributes:
        band_struc: Object which contains band structure informations.
    '''

    def __init__(self, band_struc):
        '''Creates NXY-instance.

        Args:
            band_struc: Object which contains band structure informations.
        '''
        self.band_struc = band_struc


    @classmethod
    def fromfile(cls, fobj):
        '''Creates an XYZ instance from a file.

        Args:
            fobj: filename or file like object containing band structure
            informations in NXY-format.

        '''
        band_struc = np.loadtxt(fobj, dtype=float)
        return cls(band_struc)


    def tofile(self, fobj):
        '''Writes an NXY file.

        Args:
            fobj: File name or file object where band structure should be written.
        '''
        np.savetxt(self.band_struc, fobj)


    def equals(self, other, abstolerance=_ABSTOLERANCE, reltolerance=_RELTOLERANCE):
        '''Checks whether object equals to an other one.

        Args:
            other (Nxy): Other Nxy object.
            tolerance (float): Maximal allowed deviation in floating point
                numbers (e.g. coordinates).
        '''
        if not np.allclose(self.band_struc, other.band_struc, rtol=reltolerance,
                           atol=abstolerance):
            return False
        return True
